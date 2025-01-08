import os
import csv
import json
import numpy as np
from sklearn.manifold import MDS
import plotly.graph_objects as go
import networkx as nx
import subprocess
import re
from http.server import HTTPServer, SimpleHTTPRequestHandler, ThreadingHTTPServer
import threading
import webbrowser
import urllib.parse
from tqdm import tqdm
import time
import sys

def extract_ref_number(strain_info):
    """Extract RSP number from strain info string"""
    match = re.search(r'RSP\d+', strain_info)
    return match.group(0) if match else None

def load_strain_data(folder_path):
    """Load genetic relationship data from all strain folders and their relationships"""
    strains_data = {}
    all_relationships = set()
    
    # Walk through all subdirectories
    for root, dirs, files in os.walk(folder_path):
        for dir in dirs:
            if not dir.startswith('.'):  # Skip hidden directories
                # Clean strain name by removing extra spaces
                strain_name = ' '.join(dir.split('-')[0].strip().split())
                
                # Check if this is a successfully scraped strain
                metadata_file = os.path.join(root, dir, f"{strain_name.replace(' ', '_')}.metadata.csv")
                chemicals_file = os.path.join(root, dir, f"{strain_name.replace(' ', '_')}.chemicals.csv")
                variants_file = os.path.join(root, dir, f"{strain_name.replace(' ', '_')}.variants.csv")
                
                # Extract RSP number from directory name
                rsp_match = re.search(r'-rsp(\d+)', dir.lower())
                rsp = f"RSP{rsp_match.group(1)}" if rsp_match else ''
                
                # Mark strain as complete if all files exist and have data
                is_complete = all([
                    os.path.exists(f) and os.path.getsize(f) > 0 
                    for f in [metadata_file, chemicals_file, variants_file]
                ])
                
                strains_data[strain_name] = {
                    'complete': is_complete,
                    'rsp': rsp,
                    'dir_name': dir
                }
                
                # Add relationships if they exist
                if os.path.exists(variants_file):
                    with open(variants_file, 'r') as f:
                        reader = csv.DictReader(f)
                        for row in reader:
                            if row.get('Distance') and row.get('Strain'):
                                # Clean relationship strain name
                                rel_strain = ' '.join(row['Strain'].strip().split())
                                rel = (strain_name, rel_strain, float(row['Distance']))
                                all_relationships.add(rel)
                                
                                # Add related strain to strains_data if not exists
                                if rel_strain not in strains_data:
                                    strains_data[rel_strain] = {
                                        'complete': False,
                                        'rsp': row.get('RSP', ''),
                                        'dir_name': ''
                                    }
                
                # Add terpene data if available
                if os.path.exists(chemicals_file):
                    with open(chemicals_file, 'r') as f:
                        reader = csv.DictReader(f)
                        terpenes = {}
                        for row in reader:
                            name = row.get('Name', '').lower()
                            if ('terpene' in name or 
                                any(t in name for t in ['myrcene', 'limonene', 'pinene', 'caryophyllene'])):
                                # Strip percentage sign and convert to float
                                value = row.get('Value', '0')
                                value = value.strip().rstrip('%')  # Remove % sign and whitespace
                                try:
                                    terpenes[row['Name']] = float(value)
                                except ValueError:
                                    print(f"Warning: Could not convert value '{value}' to float for terpene {row['Name']}")
                                    terpenes[row['Name']] = 0.0
                        strains_data[strain_name]['terpenes'] = terpenes
    
    return strains_data, all_relationships

def parse_summary_relationships(content):
    """Parse relationships from summary file content"""
    relationships = {
        'all_samples': {},
        'base_tree': {},
        'all_strains': set()
    }
    
    # Find the sections
    all_samples_match = re.search(r'Nearest Genetic Relatives \(All Samples\):\n(.*?)\n\n', content, re.DOTALL)
    base_tree_match = re.search(r'Nearest Genetic Relatives \(Base Tree\):\n(.*?)\n\n', content, re.DOTALL)
    
    if all_samples_match:
        for line in all_samples_match.group(1).strip().split('\n'):
            if line.strip():
                parts = line.strip().split(' - ', 1)
                if len(parts) == 2:
                    distance = float(parts[0].strip())
                    strain_info = parts[1].strip()
                    # Store the full strain info including RSP number
                    strain_name = strain_info.split('(')[0].strip()
                    rsp_number = re.search(r'\((RSP\d+)\)', strain_info)
                    if rsp_number:
                        relationships['all_samples'][strain_name] = {
                            'distance': distance,
                            'rsp': rsp_number.group(1)
                        }
                        relationships['all_strains'].add(strain_info)
    
    if base_tree_match:
        for line in base_tree_match.group(1).strip().split('\n'):
            if line.strip():
                parts = line.strip().split(' - ', 1)
                if len(parts) == 2:
                    distance = float(parts[0].strip())
                    strain_info = parts[1].strip()
                    strain_name = strain_info.split('(')[0].strip()
                    rsp_number = re.search(r'\((RSP\d+)\)', strain_info)
                    if rsp_number:
                        relationships['base_tree'][strain_name] = {
                            'distance': distance,
                            'rsp': rsp_number.group(1)
                        }
                        relationships['all_strains'].add(strain_info)
    
    return relationships

def create_distance_matrix(strains_data, all_relationships):
    """Create a distance matrix including all known relationships"""
    # Create a list of all strain names
    all_strain_names = set()
    
    # Add names from strains_data
    for strain_name in strains_data.keys():
        all_strain_names.add(strain_name)
    
    # Add names from relationships
    for strain1, strain2, _ in all_relationships:
        all_strain_names.add(strain1)
        all_strain_names.add(strain2)
    
    strain_names = list(all_strain_names)
    n = len(strain_names)
    
    # Initialize distance matrix with max distance (1.0)
    distances = np.ones((n, n))
    np.fill_diagonal(distances, 0)  # Set diagonal to 0
    
    # Fill in known distances from relationships
    name_to_index = {name: i for i, name in enumerate(strain_names)}
    for strain1, strain2, distance in all_relationships:
        i = name_to_index[strain1]
        j = name_to_index[strain2]
        distances[i,j] = distance
        distances[j,i] = distance  # Make matrix symmetric
    
    return distances, strain_names

def scrape_missing_strain(strain_info):
    """Scrape data for a missing strain"""
    ref_number = extract_ref_number(strain_info)
    if ref_number:
        ref_number = ref_number.lower()
        print(f"Scraping data for {strain_info}...")
        try:
            subprocess.run(['python', 'kaana_scraper.py', '-u', ref_number], check=True)
            print(f"Successfully scraped {strain_info}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error scraping {strain_info}: {e}")
            return False
    return False

def create_2d_visualization(strains_data, all_relationships, terpene_relationships):
    """Create interactive visualization using Vis.js"""
    # Create nodes and edges for Vis.js
    nodes = []
    relationships = []
    
    # Track RSP numbers to avoid duplicates
    seen_rsp = {}
    
    # First pass: Find completed nodes for each RSP
    for strain_name, data in strains_data.items():
        rsp = data.get('rsp', '').upper()
        if rsp:  # Only track nodes that have an RSP number
            if rsp not in seen_rsp:
                seen_rsp[rsp] = {
                    'name': strain_name,
                    'complete': data['complete']
                }
            elif data['complete'] and not seen_rsp[rsp]['complete']:
                # Replace existing incomplete node with complete one
                seen_rsp[rsp] = {
                    'name': strain_name,
                    'complete': True
                }

    # Second pass: Create nodes, skipping duplicates
    for strain_name, data in strains_data.items():
        rsp = data.get('rsp', '').upper()
        
        # Skip if this is a duplicate RSP and not the chosen representative
        if rsp and seen_rsp[rsp]['name'] != strain_name:
            continue
            
        nodes.append({
            'id': strain_name,
            'label': strain_name,
            'title': f"{strain_name}<br>RSP: {data.get('rsp', '')}<br>{'Has full data' if data['complete'] else 'Click to scrape data'}",
            'color': {
                'background': '#2B7CE9' if data['complete'] else '#cccccc',
                'border': '#2B7CE9' if data['complete'] else '#666666'
            },
            'rsp': data.get('rsp', ''),
            'complete': data['complete']
        })
    
    # Update relationships to use the chosen strain names
    for strain1, strain2, distance in all_relationships:
        # Get the representative strain names if they exist
        rsp1 = strains_data.get(strain1, {}).get('rsp', '').upper()
        rsp2 = strains_data.get(strain2, {}).get('rsp', '').upper()
        
        from_strain = seen_rsp.get(rsp1, {}).get('name', strain1) if rsp1 else strain1
        to_strain = seen_rsp.get(rsp2, {}).get('name', strain2) if rsp2 else strain2
        
        relationships.append({
            'from': from_strain,
            'to': to_strain,
            'distance': distance
        })

    # Read the HTML template
    with open('visualization_template.html', 'r', encoding='utf-8') as f:
        template = f.read()
    
    # Convert data to JSON and properly escape it
    nodes_json = json.dumps(nodes)
    relationships_json = json.dumps(relationships)
    terpene_relationships_json = json.dumps(terpene_relationships)
    
    # Create a data initialization script
    data_script = f"""
        window.INITIAL_DATA = {{
            nodes: {nodes_json},
            relationships: {relationships_json},
            terpeneRelationships: {terpene_relationships_json}
        }};
    """
    
    # Insert the data into the template
    html_content = template.replace(
        '{{DATA_INITIALIZATION}}',
        data_script
    )
    
    return html_content

def calculate_terpene_relationships(strains_data):
    """Calculate similarity relationships between strains based on their terpene profiles"""
    terpene_relationships = []
    
    # Get all strains with terpene data
    strains_with_terpenes = {
        name: data for name, data in strains_data.items() 
        if data.get('terpenes') and data['complete']
    }
    
    # Calculate cosine similarity between each pair of strains
    for strain1, data1 in strains_with_terpenes.items():
        terpenes1 = data1['terpenes']
        for strain2, data2 in strains_with_terpenes.items():
            if strain1 >= strain2:  # Skip duplicate pairs and self-comparisons
                continue
                
            terpenes2 = data2['terpenes']
            
            # Get all unique terpenes
            all_terpenes = set(terpenes1.keys()) | set(terpenes2.keys())
            
            # Calculate cosine similarity
            dot_product = sum(terpenes1.get(t, 0) * terpenes2.get(t, 0) for t in all_terpenes)
            norm1 = sum(v * v for v in terpenes1.values()) ** 0.5
            norm2 = sum(v * v for v in terpenes2.values()) ** 0.5
            
            if norm1 and norm2:  # Avoid division by zero
                similarity = dot_product / (norm1 * norm2)
                # Convert similarity to distance (0 = identical, 1 = completely different)
                distance = 1 - similarity
                
                # Only include relationships with meaningful similarity
                if distance < 0.8:  # Adjust this threshold as needed
                    terpene_relationships.append({
                        'from': strain1,
                        'to': strain2,
                        'distance': distance
                    })
    
    return terpene_relationships

class ScraperHandler(SimpleHTTPRequestHandler):
    def get_strain_data(self, strain_name, rsp):
        """Read strain data from files"""
        try:
            print(f"\n=== Reading data for {strain_name} (RSP: {rsp}) ===")
            
            # Construct directory path - fix the path to look in current directory
            dir_name = f"{strain_name.replace(' ', '_')}-{rsp.lower()}"
            base_path = os.path.join('.', 'plants', dir_name)  # Changed path
            print(f"Looking in directory: {base_path}")
            
            if not os.path.exists(base_path):
                print(f"Directory not found: {base_path}")
                return None
            
            # Construct file paths
            base_name = strain_name.replace(' ', '_')
            summary_file = os.path.join(base_path, f"{base_name}_summary.txt")
            chemicals_file = os.path.join(base_path, f"{base_name}.chemicals.csv")
            metadata_file = os.path.join(base_path, f"{base_name}.metadata.csv")
            
            print(f"Checking files exist:")
            print(f"Summary: {os.path.exists(summary_file)}")
            print(f"Chemicals: {os.path.exists(chemicals_file)}")
            print(f"Metadata: {os.path.exists(metadata_file)}")
            
            # Read summary file
            with open(summary_file, 'r', encoding='utf-8') as f:
                summary = f.read()
                print("✓ Read summary file")
            
            # Read chemicals file
            chemicals = []
            with open(chemicals_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                chemicals = [
                    {
                        'Name': row['Name'],
                        'Value': row['Value'],
                        'Type': 'Cannabinoid' if 'THC' in row['Name'] or 'CBD' in row['Name'] 
                               else 'Terpenoid'
                    }
                    for row in reader
                ]
                print(f"✓ Read {len(chemicals)} chemical entries")
            
            # Read metadata file
            with open(metadata_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                metadata = next(reader)
                print("✓ Read metadata")
            
            data = {
                'summary': summary,
                'chemicals': chemicals,
                'metadata': metadata
            }
            print("=== Successfully loaded all data ===\n")
            return data
            
        except Exception as e:
            print(f"\n!!! Error reading strain data: {str(e)}")
            import traceback
            traceback.print_exc()
            return None

    def do_GET(self):
        print(f"\nGET request: {self.path}")
        
        if self.path == '/':
            # Redirect root to visualization.html
            self.send_response(301)
            self.send_header('Location', '/visualization.html')
            self.end_headers()
            return
            
        elif self.path == '/visualization.html':
            try:
                with open('visualization.html', 'rb') as f:
                    content = f.read()
                    self.send_response(200)
                    self.send_header('Content-type', 'text/html')
                    self.end_headers()
                    self.wfile.write(content)
                    print("✓ Served visualization.html")
                    return
            except Exception as e:
                print(f"!!! Error serving visualization: {e}")
                self.send_error(500)
                return
                
        elif self.path.startswith('/scrape/'):
            try:
                rsp = self.path.split('/scrape/')[1]
                print(f"Scraping strain with RSP: {rsp}")
                
                # Get the absolute path to kaana_scraper.py
                scraper_path = os.path.join(os.path.dirname(__file__), 'kaana_scraper.py')
                
                # Run the scraper
                process = subprocess.run(
                    [sys.executable, scraper_path, '-u', rsp],
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                if process.returncode == 0:
                    # Find the newly created directory
                    strain_dir = None
                    strain_name = None
                    
                    for root, dirs, files in os.walk('./plants'):
                        for dir in dirs:
                            if rsp.lower() in dir.lower():
                                strain_dir = dir
                                strain_name = ' '.join(dir.split('-')[0].strip().split('_'))
                                break
                    
                    if strain_dir and strain_name:
                        # Get the strain data
                        strain_data = self.get_strain_data(strain_name, rsp)
                        
                        self.send_response(200)
                        self.send_header('Content-type', 'application/json')
                        self.send_header('Access-Control-Allow-Origin', '*')
                        self.end_headers()
                        self.wfile.write(json.dumps({
                            'success': True,
                            'strain_name': strain_name,
                            'strain_data': strain_data
                        }).encode())
                    else:
                        raise Exception("Could not find scraped strain directory")
                else:
                    raise Exception(f"Scraper failed: {process.stderr}")
                    
            except Exception as e:
                print(f"Error during scraping: {str(e)}")
                self.send_response(500)
                self.send_header('Content-type', 'application/json')
                self.send_header('Access-Control-Allow-Origin', '*')
                self.end_headers()
                self.wfile.write(json.dumps({
                    'success': False,
                    'error': str(e)
                }).encode())
                return
                
        elif self.path.startswith('/strain_data/'):
            try:
                path_parts = self.path.split('/strain_data/')[1]
                strain_name, rsp = urllib.parse.unquote(path_parts).split('|')
                print(f"\nRequested data for: {strain_name} (RSP: {rsp})")
                
                strain_data = self.get_strain_data(strain_name, rsp)
                if strain_data:
                    self.send_response(200)
                    self.send_header('Content-type', 'application/json')
                    self.send_header('Access-Control-Allow-Origin', '*')
                    self.end_headers()
                    self.wfile.write(json.dumps({
                        'success': True,
                        'data': strain_data
                    }).encode())
                    print("✓ Sent strain data")
                else:
                    raise Exception("Could not read strain data")
                
            except Exception as e:
                print(f"!!! Error: {str(e)}")
                self.send_response(500)
                self.send_header('Content-type', 'application/json')
                self.send_header('Access-Control-Allow-Origin', '*')
                self.end_headers()
                self.wfile.write(json.dumps({
                    'success': False,
                    'error': str(e)
                }).encode())
                
        else:
            # Handle other routes as before
            return super().do_GET()

def start_server(port=8000):
    """Start the HTTP server"""
    server = HTTPServer(('localhost', port), ScraperHandler)
    server_thread = threading.Thread(target=server.serve_forever)
    server_thread.daemon = True
    server_thread.start()
    return server

def main():
    print("\n=== Starting Visualization Server ===")
    print("Loading strain data...")
    strains_data, all_relationships = load_strain_data('.')
    print(f"Loaded data for {len(strains_data)} strains")
    print(f"Found {len(all_relationships)} total relationships")
    
    print("\nCalculating terpene relationships...")
    terpene_relationships = calculate_terpene_relationships(strains_data)
    print(f"Found {len(terpene_relationships)} terpene relationships")
    
    print("\nCreating visualization...")
    html_content = create_2d_visualization(strains_data, all_relationships, terpene_relationships)
    
    # Save the HTML file
    print("\nSaving HTML file...")
    with open('visualization.html', 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    # Start the server
    port = 8000
    server = ThreadingHTTPServer(('', port), ScraperHandler)
    print(f"\nServer started at http://localhost:{port}")
    print("Press Ctrl+C to stop the server")
    
    # Open the visualization in the default browser
    webbrowser.open(f'http://localhost:{port}/visualization.html')
    
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down server...")
        server.shutdown()
        server.server_close()

if __name__ == "__main__":
    main() 