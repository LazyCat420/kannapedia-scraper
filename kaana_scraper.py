import requests
import os
import argparse
from bs4 import BeautifulSoup
from playwright.async_api import async_playwright
import asyncio
import re
import csv

async def scrape_strain_data(rsp_number):
    """Scrape data for a specific strain using its RSP number"""
    print(f"Starting scrape for {rsp_number}")
    base_url = "https://www.kannapedia.net/strains/"
    url = f"{base_url}{rsp_number}"
    
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        try:
            print(f"Loading page: {url}")
            await page.goto(url, wait_until='networkidle')
            await page.wait_for_selector('h1.StrainInfo--title', timeout=30000)
            
            # Extract all data using JavaScript evaluation
            strain_data = await page.evaluate("""
            () => {
                    const data = {
                        name: '',
                        general_info: {},
                        chemical_content: {
                            cannabinoids: {},
                            terpenoids: {}
                        },
                        genetic_relationships: {
                            all_samples: [],
                            base_tree: [],
                            most_distant: []
                        },
                        blockchain: {}
                    };
                    
                    // Get strain name
                    const titleElem = document.querySelector('h1.StrainInfo--title');
                    data.name = titleElem ? titleElem.textContent.trim() : '';
                    console.log('Found strain name:', data.name);
                    
                    // Get general information
                    const generalSection = Array.from(document.querySelectorAll('h2')).find(
                        h2 => h2.textContent.includes('General Information')
                    )?.parentElement;
                    
                    if (generalSection) {
                        const rows = generalSection.querySelectorAll('dt, dd');
                        for (let i = 0; i < rows.length; i += 2) {
                            if (rows[i] && rows[i + 1]) {
                                const key = rows[i].textContent.trim();
                                const value = rows[i + 1].textContent.trim();
                                data.general_info[key] = value;
                            }
                        }
                    }
                    
                    // Get grower information
                    const growerText = document.querySelector('.StrainInfo--grower');
                    if (growerText) {
                        data.general_info['Grower'] = growerText.textContent.replace('Grower:', '').trim();
                    }
                    
                    // Get chemical content
                    const chemicalSection = Array.from(document.querySelectorAll('h2')).find(
                        h2 => h2.textContent.includes('Chemical Information')
                    )?.parentElement;
                    
                    if (chemicalSection) {
                        // Get Cannabinoids
                        const cannabinoidSection = Array.from(chemicalSection.querySelectorAll('h3')).find(
                            h3 => h3.textContent.includes('Cannabinoids')
                        );
                        if (cannabinoidSection) {
                            const items = Array.from(cannabinoidSection.parentElement.querySelectorAll('dt, dd'));
                            for (let i = 0; i < items.length; i += 2) {
                                if (items[i] && items[i + 1]) {
                                    const name = items[i].textContent.trim();
                                    const value = items[i + 1].textContent.trim();
                                    if (value !== 'n/a' && !value.toLowerCase().includes('no information')) {
                                        data.chemical_content.cannabinoids[name] = value;
                                    }
                                }
                            }
                        }
                        
                        // Get Terpenoids
                        const terpenoidSection = Array.from(chemicalSection.querySelectorAll('h3')).find(
                            h3 => h3.textContent.includes('Terpenoids')
                        );
                        if (terpenoidSection) {
                            const items = Array.from(terpenoidSection.parentElement.querySelectorAll('dt, dd'));
                            for (let i = 0; i < items.length; i += 2) {
                                if (items[i] && items[i + 1]) {
                                    const name = items[i].textContent.trim();
                                    const value = items[i + 1].textContent.trim();
                                    if (value !== 'n/a' && !value.toLowerCase().includes('no information')) {
                                        data.chemical_content.terpenoids[name] = value;
                                    }
                                }
                            }
                        }
                    }
                    
                    // Get heterozygosity
                    const heteroText = document.evaluate(
                        "//text()[contains(., 'Heterozygosity:')]",
                    document,
                    null,
                    XPathResult.FIRST_ORDERED_NODE_TYPE,
                    null
                    ).singleNodeValue;
                    if (heteroText) {
                        const match = heteroText.textContent.match(/Heterozygosity:\s*([\d.]+%)/);
                        if (match) {
                            data.general_info['Reported Heterozygosity'] = match[1];
                        }
                    }
                    
                    // Get rarity
                    const rarityText = document.evaluate(
                        "//text()[contains(., 'Rarity:')]",
                    document,
                    null,
                    XPathResult.FIRST_ORDERED_NODE_TYPE,
                    null
                    ).singleNodeValue;
                    if (rarityText) {
                        const match = rarityText.textContent.match(/Rarity:\s*(\w+)/);
                        if (match) {
                            data.general_info['Rarity'] = match[1];
                        }
                    }
                    
                    // Extract genetic relationships
                    const extractRelationships = (title) => {
                        const results = [];
                        const listItems = Array.from(document.querySelectorAll('li')).filter(li => {
                            const text = li.textContent.trim();
                            return text.match(/^\d+\.\d+\s+.+\(RSP\d+\)/);
                        });
                        
                        listItems.forEach(li => {
                            const text = li.textContent.trim();
                            const match = text.match(/^(\d+\.\d+)\s+(.+?)\s*\((RSP\d+)\)/i);
                            if (match) {
                                results.push({
                                    distance: parseFloat(match[1]),
                                    strain: match[2].trim(),
                                    rsp: match[3].toLowerCase()
                                });
                            }
                        });
                        return results;
                    };
                    
                    // Get all genetic relationships
                    const geneticSections = document.querySelectorAll('h3');
                    geneticSections.forEach(section => {
                        const title = section.textContent.trim().toLowerCase();
                        if (title.includes('all samples')) {
                            data.genetic_relationships.all_samples = extractRelationships(title);
                        } else if (title.includes('base tree')) {
                            data.genetic_relationships.base_tree = extractRelationships(title);
                        } else if (title.includes('most genetically distant')) {
                            data.genetic_relationships.most_distant = extractRelationships(title);
                        }
                    });
                    console.log('Found genetic relationships:', 
                        Object.keys(data.genetic_relationships).map(k => 
                            `${k}: ${data.genetic_relationships[k].length} items`
                        )
                    );
                    
                    // Get blockchain information
                    const txidElem = document.evaluate(
                        "//dt[contains(text(), 'Transaction ID')]/following-sibling::dd[1]",
                        document,
                        null,
                        XPathResult.FIRST_ORDERED_NODE_TYPE,
                        null
                    ).singleNodeValue;
                    if (txidElem) {
                        data.blockchain.txid = txidElem.textContent.trim();
                    }
                    
                    const shasumElem = document.evaluate(
                        "//dt[contains(text(), 'SHASUM Hash')]/following-sibling::dd[1]",
                        document,
                        null,
                        XPathResult.FIRST_ORDERED_NODE_TYPE,
                        null
                    ).singleNodeValue;
                    if (shasumElem) {
                        data.blockchain.shasum = shasumElem.textContent.trim();
                    }
                    console.log('Found blockchain info:', data.blockchain);
                    
                    console.log('Found general info:', data.general_info);
                    console.log('Found chemical content:', data.chemical_content);
                    
                    return data;
                }
            """)
            
            print("Extracted data:", strain_data)
            
            # Create directory structure
            strain_dir = os.path.join('plants', f"{strain_data['name'].replace(' ', '_')}-{rsp_number}")
            os.makedirs(strain_dir, exist_ok=True)
            
            # Save metadata CSV
            metadata_path = os.path.join(strain_dir, f"{strain_data['name'].replace(' ', '_')}.metadata.csv")
            with open(metadata_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(['Field', 'Value'])
                for key, value in strain_data['general_info'].items():
                    writer.writerow([key, value])
            
            # Save chemicals CSV
            chemicals_path = os.path.join(strain_dir, f"{strain_data['name'].replace(' ', '_')}.chemicals.csv")
            with open(chemicals_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(['Type', 'Name', 'Value'])
                for name, value in strain_data['chemical_content']['cannabinoids'].items():
                    writer.writerow(['Cannabinoid', name, value])
                for name, value in strain_data['chemical_content']['terpenoids'].items():
                    writer.writerow(['Terpenoid', name, value])
            
            # Save variants CSV
            variants_path = os.path.join(strain_dir, f"{strain_data['name'].replace(' ', '_')}.variants.csv")
            with open(variants_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(['Type', 'Distance', 'Strain', 'RSP'])
                for rel_type in ['all_samples', 'base_tree', 'most_distant']:
                    for rel in strain_data['genetic_relationships'][rel_type]:
                        writer.writerow([rel_type, rel['distance'], rel['strain'], rel['rsp'].upper()])
            
            # Save summary text file
            summary_path = os.path.join(strain_dir, f"{strain_data['name'].replace(' ', '_')}_summary.txt")
            with open(summary_path, 'w', encoding='utf-8') as f:
        f.write(f"{'='*80}\n")
                f.write(f"{strain_data['name']} ({rsp_number.upper()}) Summary\n")
        f.write(f"{'='*80}\n\n")
        
                # Write general information
        f.write("GENERAL INFORMATION\n")
        f.write(f"{'-'*80}\n")
                for key, value in strain_data['general_info'].items():
                    f.write(f"{key}: {value}\n")
                f.write("\n")
                
                # Write chemical content
        f.write("CHEMICAL CONTENT\n")
        f.write(f"{'-'*80}\n")
        f.write("Cannabinoids:\n")
                for name, value in strain_data['chemical_content']['cannabinoids'].items():
                    f.write(f"  {name}: {value}\n")
        f.write("\nTerpenoids:\n")
                for name, value in strain_data['chemical_content']['terpenoids'].items():
                    f.write(f"  {name}: {value}\n")
        f.write("\n")
        
                # Write genetic relationships
        f.write("GENETIC RELATIONSHIPS\n")
        f.write(f"{'-'*80}\n")
        
        f.write("Nearest Genetic Relatives (All Samples):\n")
                for rel in strain_data['genetic_relationships']['all_samples']:
                    f.write(f"  {rel['distance']:.3f} - {rel['strain']} ({rel['rsp'].upper()})({rel['rsp']})\n")
        f.write("\n")
        
        f.write("Nearest Genetic Relatives (Base Tree):\n")
                for rel in strain_data['genetic_relationships']['base_tree']:
                    f.write(f"  {rel['distance']:.3f} - {rel['strain']} ({rel['rsp'].upper()})({rel['rsp']})\n")
        f.write("\n")
        
                if strain_data['genetic_relationships']['most_distant']:
        f.write("Most Genetically Distant Strains:\n")
                    for rel in strain_data['genetic_relationships']['most_distant']:
                        f.write(f"  {rel['distance']:.3f} - {rel['strain']} ({rel['rsp'].upper()})({rel['rsp']})\n")
        f.write("\n")
        
                # Write blockchain information
        f.write("BLOCKCHAIN INFORMATION\n")
        f.write(f"{'-'*80}\n")
                if strain_data['blockchain'].get('txid'):
                    f.write(f"Transaction ID: {strain_data['blockchain']['txid']}\n")
                if strain_data['blockchain'].get('shasum'):
                    f.write(f"SHASUM Hash: {strain_data['blockchain']['shasum']}\n")
            
            print(f"Saved all data for {strain_data['name']}")
            await browser.close()
            return True
            
        except Exception as e:
            print(f"Error during scraping: {str(e)}")
            print("Full error:")
            import traceback
            traceback.print_exc()
            await browser.close()
            raise

def main():
    parser = argparse.ArgumentParser(description="Scrape strain data from Kannapedia")
    parser.add_argument('-u', '--url', required=True, help='RSP number of the strain to scrape')
    args = parser.parse_args()
    
    # Clean the RSP number
    rsp_number = args.url.lower()
    if not rsp_number.startswith('rsp'):
        rsp_number = 'rsp' + rsp_number.replace('rsp', '')
    
    try:
        asyncio.run(scrape_strain_data(rsp_number))
        print("Scraping completed successfully")
    except Exception as e:
        print(f"Error during scraping: {str(e)}")
        raise

if __name__ == "__main__":
    main()