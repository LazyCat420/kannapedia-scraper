import os
import csv
import json
import numpy as np
from sklearn.manifold import MDS
import plotly.graph_objects as go
import networkx as nx
import subprocess
import re
from http.server import HTTPServer, SimpleHTTPRequestHandler
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
                # Find the summary txt file
                strain_name = dir.split('-')[0].strip()
                summary_file = f"{strain_name}_summary.txt"
                file_path = os.path.join(folder_path, dir, summary_file)
                
                if os.path.exists(file_path):
                    with open(file_path, 'r') as f:
                        content = f.read()
                        
                        # Extract RSP number
                        rsp_match = re.search(r'\((RSP\d+)\)', content)
                        ref_number = rsp_match.group(1) if rsp_match else None
                        
                        # Extract relationships
                        relationships = parse_summary_relationships(content)
                        
                        strains_data[strain_name] = {
                            'all_samples': relationships['all_samples'],
                            'base_tree': relationships['base_tree'],
                            'ref_number': ref_number
                        }
                        
                        # Add all related strains to the set
                        for strain_info in relationships['all_strains']:
                            all_relationships.add(strain_info)
    
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

def create_distance_matrix(strains_data, all_strains):
    """Create a distance matrix including all known relationships"""
    # Create a list of all strain names
    all_strain_names = set()
    for strain_info in all_strains:
        strain_name = strain_info.split('(')[0].strip()
        all_strain_names.add(strain_name)
    for strain_name in strains_data.keys():
        all_strain_names.add(strain_name)
    
    strain_names = list(all_strain_names)
    n = len(strain_names)
    
    # Initialize distance matrix with max distance (1.0)
    distances = np.ones((n, n))
    np.fill_diagonal(distances, 0)  # Set diagonal to 0
    
    # Fill in known distances
    for i, strain1 in enumerate(strain_names):
        for j, strain2 in enumerate(strain_names):
            if strain1 in strains_data and strain2 in strains_data[strain1]['all_samples']:
                # Extract the distance from the dictionary structure
                distances[i,j] = strains_data[strain1]['all_samples'][strain2]['distance']
                distances[j,i] = distances[i,j]  # Make matrix symmetric
    
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

def create_3d_visualization(strains_data, all_strains):
    """Create 3D visualization using MDS"""
    # Create distance matrix
    distances, strain_names = create_distance_matrix(strains_data, all_strains)
    
    # Apply MDS to get 3D coordinates
    mds = MDS(n_components=3, dissimilarity='precomputed', random_state=42)
    coords = mds.fit_transform(distances)
    
    # Create graph for edge visualization
    G = nx.Graph()
    
    # Add edges for closely related strains (distance < 0.15)
    for i, strain1 in enumerate(strain_names):
        for j, strain2 in enumerate(strain_names):
            if distances[i,j] < 0.15 and i != j:
                G.add_edge(strain1, strain2, weight=distances[i,j])
    
    # Create edge traces
    edge_x = []
    edge_y = []
    edge_z = []
    
    for edge in G.edges():
        x0, y0, z0 = coords[strain_names.index(edge[0])]
        x1, y1, z1 = coords[strain_names.index(edge[1])]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        edge_z.extend([z0, z1, None])
    
    # Create customdata with RSP numbers
    customdata = []
    for name in strain_names:
        if name not in strains_data:
            # Look for RSP number in relationships
            rsp_info = None
            # First check in all_samples relationships
            for strain1 in strains_data:
                if name in strains_data[strain1]['all_samples']:
                    rsp_info = strains_data[strain1]['all_samples'][name]['rsp']
                    break
            # If not found, check in base_tree relationships
            if not rsp_info:
                for strain1 in strains_data:
                    if name in strains_data[strain1]['base_tree']:
                        rsp_info = strains_data[strain1]['base_tree'][name]['rsp']
                        break
            customdata.append(rsp_info)
        else:
            # For strains we have data for, use their ref_number
            customdata.append(strains_data[name]['ref_number'])
    
    # Determine node colors based on whether we have full data
    node_colors = ['darkblue' if name in strains_data else 'lightgray' for name in strain_names]
    
    # Create visualization
    fig = go.Figure(data=[
        # Add edges
        go.Scatter3d(
            x=edge_x, y=edge_y, z=edge_z,
            mode='lines',
            line=dict(color='lightgray', width=1),
            hoverinfo='none'
        ),
        # Add nodes
        go.Scatter3d(
            x=coords[:,0],
            y=coords[:,1],
            z=coords[:,2],
            mode='markers+text',
            text=strain_names,
            hovertext=[f"{name}<br>RSP: {customdata[i] if customdata[i] else 'N/A'}<br>{'Has full data' if name in strains_data else 'Click to scrape data'}" 
                      for i, name in enumerate(strain_names)],
            hoverinfo='text',
            marker=dict(
                size=8,
                color=node_colors,
                opacity=0.8
            ),
            customdata=customdata
        )
    ])
    
    # Update layout
    fig.update_layout(
        title='3D Cannabis Strain Genetic Relationships',
        scene=dict(
            xaxis_title='Dimension 1',
            yaxis_title='Dimension 2',
            zaxis_title='Dimension 3'
        ),
        showlegend=False,
        margin=dict(l=0, r=0, t=30, b=0)
    )
    
    return fig

class ScraperHandler(SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path.startswith('/scrape/'):
            print("\n" + "="*50)
            print("SCRAPING REQUEST STARTED")
            print("="*50)
            
            try:
                # Get the Python executable from the current virtual environment
                python_executable = sys.executable
                print(f"Using Python from: {python_executable}")

                # Extract RSP number from path
                rsp_number = self.path.split('/scrape/')[1]
                rsp_number = urllib.parse.unquote(rsp_number)
                print(f"[1] Attempting to scrape RSP number: {rsp_number}")
                
                # Run the scraper script using the virtual environment's Python
                cmd = [python_executable, 'kaana_scraper.py', '-u', rsp_number]
                print(f"[2] Running command: {' '.join(cmd)}")
                
                # First ensure playwright is installed
                try:
                    import playwright
                except ImportError:
                    print("[2a] Installing playwright...")
                    subprocess.run([python_executable, '-m', 'pip', 'install', 'playwright'], check=True)
                    print("[2b] Installing playwright browsers...")
                    subprocess.run([python_executable, '-m', 'playwright', 'install'], check=True)
                
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
                
                stdout, stderr = process.communicate()
                
                print("[3] Scraper Output:")
                print(stdout)
                if stderr:
                    print("[3a] Scraper Errors:")
                    print(stderr)
                
                if process.returncode == 0:
                    print("[4] Scraping completed successfully")
                    self.send_response(200)
                    self.send_header('Content-type', 'application/json')
                    self.send_header('Access-Control-Allow-Origin', '*')
                    self.end_headers()
                    self.wfile.write(json.dumps({
                        'success': True,
                        'message': 'Scraping completed successfully',
                        'output': stdout
                    }).encode())
                else:
                    raise Exception(f"Scraper failed with return code {process.returncode}")
                
            except Exception as e:
                print(f"[ERROR] Scraping failed: {str(e)}")
                self.send_response(500)
                self.send_header('Content-type', 'application/json')
                self.send_header('Access-Control-Allow-Origin', '*')
                self.end_headers()
                self.wfile.write(json.dumps({
                    'success': False,
                    'error': str(e)
                }).encode())
            
            print("="*50)
            print("SCRAPING REQUEST ENDED")
            print("="*50)
            
        else:
            return SimpleHTTPRequestHandler.do_GET(self)

def start_server(port=8000):
    """Start the HTTP server"""
    server = HTTPServer(('localhost', port), ScraperHandler)
    server_thread = threading.Thread(target=server.serve_forever)
    server_thread.daemon = True
    server_thread.start()
    return server

def main():
    # Load data and create visualization
    print("\n=== Starting Visualization Server ===")
    print("Loading strain data...")
    strains_data, all_relationships = load_strain_data('.')
    print(f"Loaded data for {len(strains_data)} strains")
    print(f"Found {len(all_relationships)} total relationships")
    
    print("\nCreating visualization...")
    fig = create_3d_visualization(strains_data, all_relationships)
    
    print("\nSaving HTML file...")
    html_content = fig.to_html(
        include_plotlyjs=True,
        full_html=True,
        include_mathjax=False,
    )
    
    # Add custom JavaScript for click handling
    click_script = """
    <script>
        function showStatus(message, isError = false) {
            const status = document.getElementById('scrape-status') || document.createElement('div');
            status.id = 'scrape-status';
            status.style.position = 'fixed';
            status.style.top = '10px';
            status.style.left = '50%';
            status.style.transform = 'translateX(-50%)';
            status.style.padding = '15px';
            status.style.backgroundColor = isError ? '#ffebee' : '#fff';
            status.style.border = `1px solid ${isError ? '#ef5350' : '#ccc'}`;
            status.style.borderRadius = '5px';
            status.style.zIndex = '1000';
            status.innerHTML = message;
            
            if (!status.parentElement) {
                document.body.appendChild(status);
            }
            
            return status;
        }

        var plot = document.getElementsByClassName('plotly-graph-div')[0];
        plot.on('plotly_click', function(data) {
            var point = data.points[0];
            if (point.customdata) {
                var rspNumber = point.customdata;
                console.log('Checking RSP number:', rspNumber);
                
                if (rspNumber) {
                    if (confirm('Would you like to scrape data for ' + point.text + ' (' + rspNumber + ')?')) {
                        const status = showStatus('Scraping data for ' + point.text + '...');
                        console.log('Starting scrape request for:', rspNumber);
                        
                        fetch('/scrape/' + rspNumber.toLowerCase())
                            .then(response => {
                                console.log('Received response:', response);
                                return response.json();
                            })
                            .then(data => {
                                console.log('Scraping result:', data);
                                if (data.success) {
                                    status.innerHTML = 'Successfully scraped data! Refreshing...';
                                    setTimeout(() => {
                                        location.reload();
                                    }, 2000);
                                } else {
                                    showStatus('Error: ' + data.error, true);
                                    setTimeout(() => {
                                        status.remove();
                                    }, 5000);
                                }
                            })
                            .catch(error => {
                                console.error('Fetch error:', error);
                                showStatus('Error: ' + error, true);
                                setTimeout(() => {
                                    status.remove();
                                }, 5000);
                            });
                    }
                } else {
                    console.error('No RSP number found for:', point.text);
                    showStatus('Error: Could not find RSP number for this strain', true);
                }
            }
        });
    </script>
    """
    
    html_content = html_content.replace('</body>', f'{click_script}</body>')
    
    with open('genetic_relationships_3d.html', 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print("\nStarting local server...")
    server = start_server()
    
    print("\nOpening visualization in browser...")
    webbrowser.open('http://localhost:8000/genetic_relationships_3d.html')
    
    print("\n=== Server Ready ===")
    print("Gray nodes indicate strains with missing data - click them to scrape")
    print("\nPress Ctrl+C to stop the server and exit")
    
    try:
        while True:
            input()
    except KeyboardInterrupt:
        print("\nShutting down server...")
        server.shutdown()

if __name__ == "__main__":
    main() 