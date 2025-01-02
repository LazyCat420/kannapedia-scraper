import requests, os, csv, argparse, asyncio
from bs4 import BeautifulSoup
from playwright.async_api import async_playwright
import re

### helpers


def download(url: str, dest_folder: str, force_redownload: bool = False):
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)  # create folder if it does not exist

    filename = url.split("/")[-1].replace(" ", "_")  # be careful with file names
    file_path = os.path.join(dest_folder, filename)

    if not force_redownload and os.path.exists(file_path):
        print(f"File {filename} already exists in folder {dest_folder}")
        return
    elif os.path.exists(file_path):
        print(f"Overwriting {file_path}")

    r = requests.get(url, stream=True)
    if r.ok:
        print("saving to", os.path.abspath(file_path))
        with open(file_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 8):
                if chunk:
                    f.write(chunk)
                    f.flush()
                    os.fsync(f.fileno())
    else:  # HTTP status code 4XX/5XX
        print("Download failed: status code {}\n{}".format(r.status_code, r.text))


### setting up

# getting arguments

parser = argparse.ArgumentParser(
    description="A program to scrape individual Kannapedia entries"
)

parser.add_argument(
    "-u",
    "--url",
    help="The strain ID to scrape (e.g. rsp11347). Required.",
    required=True
)
parser.add_argument(
    "-d",
    "--download_files",
    help="Downloads the scraped files rather than simply providing download links",
    action="store_true",
)
parser.add_argument(
    "-f",
    "--force_redownload",
    help="When downloading files, forces redownload of files that already exist",
    action="store_true",
)

args = parser.parse_args()
specific_page = args.url.strip()
download_files = args.download_files
force_redownload = args.force_redownload

# Validate the URL format
if not specific_page.startswith('rsp'):
    print("Error: URL should be a strain ID starting with 'rsp'")
    exit(1)

BASE_URL = "https://www.kannapedia.net/strains/"
HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
}

# Requesting page with error handling
print(f"Getting page for strain {specific_page}...")
try:
    page = requests.get(BASE_URL + specific_page, headers=HEADERS, timeout=10)
    page.raise_for_status()  # Raise an exception for bad status codes
except requests.RequestException as e:
    print(f"Error fetching page: {e}")
    exit(1)

soup = BeautifulSoup(page.content, "html.parser")

# Verify we got a valid page
if not soup.find(class_="StrainInfo--title"):
    print("Error: Could not find strain information on page. Please verify the strain ID.")
    exit(1)

### establishing CSV fields

# making general metadata CSV fields

metadata_fields = [
    "NAME",
    "REF NUMBER",
    "GROWER",
    "SAMPLE NAME",
    "ACCESSION DATE",
    "REPORTED SEX",
    "REPORT TYPE",
    "RARITY",
    "PLANT TYPE",
    "REPORTED HETEROZYGOSITY",
    "Y RATIO DISTRIBUTION",
    "GENETIC_RELATIVES_ALL_SAMPLES",
    "GENETIC_RELATIVES_BASE_TREE",
    "CLOSELY_RELATED_STRAINS_WITH_DISTANCE",
    "MODERATELY_RELATED_STRAINS_WITH_DISTANCE",
    "DISTANTLY_RELATED_STRAINS_WITH_DISTANCE",
    "TRANSACTION_ID",
    "SHASUM_HASH"
]

if not download_files:
    metadata_fields.append("FILES")

# making cannabinoids/terpinoids CSV fields

chemical_content_fields = [
    "THC + THCA",
    "CBD + CBDA",
    "THCV + THCVA",
    "CBC + CBCA",
    "CBG + CBGA",
    "CBN + CBNA",
    "ALPHA-BISABOLOL",
    "BORNEOL",
    "CAMPHENE",
    "CARENE",
    "CARYOPHYLLENE OXIDE",
    "BETA-CARYOPHYLLENE",
    "FENCHOL",
    "GERANIOL",
    "ALPHA-HUMULENE",
    "LIMONENE",
    "LINALOOL",
    "MYRCENE",
    "ALPHA-PHELLANDRENE",
    "TERPINOLENE",
    "ALPHA-TERPINEOL",
    "ALPHA-TERPINENE",
    "GAMMA-TERPINENE",
    "TOTAL NEROLIDOL",
    "TOTAL OCIMENE",
    "ALPHA-PINENE",
    "BETA-PINENE",
]

# making variant CSV fields

variants_fields = [
    "GENE",
    "HGVS_C",
    "HGVS_P",
    "ANNOTATION",
    "ANNOTATION IMPACT",
    "CONTIG",
    "CONTIG POS",
    "REF/ALT",
    "VAR FREQ NGS",
    "VAR FREQ C90",
]

# extracting metadata
metadata = []

strain_name_h1 = soup.find_all(class_="StrainInfo--title")[0]
strain_name = strain_name_h1.text.strip()
strain_name = (
    strain_name.replace("/", "")
    .replace("\\", "")
    .replace("(", "")
    .replace(")", "")
    .replace(".", "")
    .replace("#", "")
    .replace("'", "")
)
metadata.append(strain_name)

ref_number_p = soup.find(class_="StrainInfo--reportId")
ref_number = "".join(ref_number_p.text.split())
metadata.append(ref_number)

registrant_p = soup.find(class_="StrainInfo--registrant")
grower = registrant_p.find("a").text.strip()
metadata.append(grower)

sample_name_dt = soup.find("dt", string="Sample Name")
if sample_name_dt and sample_name_dt.find_next_sibling():
    metadata.append(sample_name_dt.find_next_sibling().text.strip())
else:
    metadata.append("")

accession_dt = soup.find("dt", string="Accession Date")
if accession_dt and accession_dt.find_next_sibling():
    metadata.append(accession_dt.find_next_sibling().text.strip())
else:
    metadata.append("")

reported_sex_dt = soup.find("dt", string="Reported Plant Sex")
if reported_sex_dt and reported_sex_dt.find_next_sibling():
    metadata.append(reported_sex_dt.find_next_sibling().text.strip())
else:
    metadata.append("")

report_type_dt = soup.find("dt", string="Report Type")
if report_type_dt and report_type_dt.find_next_sibling():
    metadata.append(report_type_dt.find_next_sibling().text.strip())
else:
    metadata.append("")

rarity = ""
rarity_plot = soup.find(class_="DataPlot Rarity")
if rarity_plot and rarity_plot.find("a"):
    rarity = rarity_plot.find("a").text.strip()
metadata.append(rarity)

plant_type = ""
plant_type_container = soup.find(class_="StrainGeneticInfo--basic")
if plant_type_container and plant_type_container.find("a"):
    plant_type = plant_type_container.find("a").text.strip()
metadata.append(plant_type)

reported_heterozygosity = ""
reported_heterozygosity_container = soup.find(class_="DataPlot Heterozygosity")
if reported_heterozygosity_container and reported_heterozygosity_container.find(
    "strong"
):
    reported_heterozygosity = reported_heterozygosity_container.find(
        "strong"
    ).text.strip()
metadata.append(reported_heterozygosity)

y_dist = ""
y_dist_container = soup.find(class_="DataPlot YRatio")
if y_dist_container and y_dist_container.find("strong"):
    y_dist = y_dist_container.find("strong").text.strip()
metadata.append(y_dist)

# Add this new function to get the genetic relationships data
async def get_page_data(url):
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        print("Loading page and waiting for data...")
        await page.goto(url)
        await page.wait_for_selector(".StrainGeneticInfo--distanceList")
        
        # Get All Samples relationships
        all_samples_data = await page.evaluate('''
            () => {
                const relationships = [];
                const allSamplesSection = document.evaluate(
                    "//h3[contains(text(), 'All Samples')]",
                    document,
                    null,
                    XPathResult.FIRST_ORDERED_NODE_TYPE,
                    null
                ).singleNodeValue?.closest('.StrainGeneticInfo--distanceList');
                
                if (allSamplesSection) {
                    const items = allSamplesSection.querySelectorAll('li');
                    items.forEach(item => {
                        const text = item.textContent.trim();
                        const match = text.match(/^(0\.\d+)/);
                        if (match) {
                            const distance = match[1];
                            const link = item.querySelector('a');
                            if (link) {
                                const name = link.textContent.trim();
                                const id = link.href.split('/').pop();
                                relationships.push(`${distance}:${name}(${id})`);
                            }
                        }
                    });
                }
                return relationships;
            }
        ''')
        
        # Get Base Tree relationships
        base_tree_data = await page.evaluate('''
            () => {
                const relationships = [];
                const baseTreeSection = document.evaluate(
                    "//h3[contains(text(), 'Base Tree')]",
                    document,
                    null,
                    XPathResult.FIRST_ORDERED_NODE_TYPE,
                    null
                ).singleNodeValue?.closest('.StrainGeneticInfo--distanceList');
                
                if (baseTreeSection) {
                    const items = baseTreeSection.querySelectorAll('li');
                    items.forEach(item => {
                        const text = item.textContent.trim();
                        const match = text.match(/^(0\.\d+)/);
                        if (match) {
                            const distance = match[1];
                            const link = item.querySelector('a');
                            if (link) {
                                const name = link.textContent.trim();
                                const id = link.href.split('/').pop();
                                relationships.push(`${distance}:${name}(${id})`);
                            }
                        }
                    });
                }
                return relationships;
            }
        ''')
        
        await browser.close()
        
        # Debug print
        print("\nAll Samples Relationships found:")
        for rel in all_samples_data:
            print(rel)
        
        print("\nBase Tree Relationships found:")
        for rel in base_tree_data:
            print(rel)
        
        return all_samples_data, base_tree_data

# Replace the genetic relationships extraction with:
print("Getting genetic relationships data...")
all_samples_data, base_tree_data = asyncio.run(
    get_page_data(f"{BASE_URL}{specific_page}")
)

# Add the relationships to metadata with proper formatting
metadata.append("|".join(all_samples_data))
metadata.append("|".join(base_tree_data))

# For distant strains (if needed)
metadata.append("")  # Placeholder for distant strains

# Update the blockchain information extraction
transaction_id = soup.find(string=lambda text: "Transaction ID" in str(text) if text else False)
if transaction_id:
    # Get the next sibling that contains the actual ID
    id_element = transaction_id.find_next()
    if id_element:
        metadata.append(id_element.text.strip())
else:
    metadata.append("")

shasum_hash = soup.find(string=lambda text: "SHASUM Hash" in str(text) if text else False)
if shasum_hash:
    # Get the next sibling that contains the actual hash
    hash_element = shasum_hash.find_next()
    if hash_element:
        metadata.append(hash_element.text.strip())
else:
    metadata.append("")

metadata_rows = [metadata]

# making folder

parent_directory = os.getcwd()

folder_name = (
    "_".join(strain_name.split()) + "-" + specific_page
)  # replacing spaces with _
path = os.path.join(parent_directory, folder_name)

print("Creating directory...")
if not os.path.exists(path):
    os.mkdir(path)
    print("Directory '% s' created" % folder_name)
else:
    print("Directory '% s' already exists" % folder_name)

# getting file downloads

files_a_tags = soup.find_all(class_="DownloadLink")
links = []
for f in files_a_tags:
    links.append(f["href"])
for x in range(len(links)):
    if links[x][0:4] != "http":
        links[x] = "https://www.kannapedia.net" + links[x]

if download_files:
    print("Skipping file downloads as they are not accessible...")
else:
    print("Saving download links to metafile...")
    metadata_rows[0].append(links)

# extracting cannabinoid/terpinoid data

chemical_content_rows = [[]]

# Extract cannabinoid data
cannabinoid_info = soup.find(class_="StrainChemicalInfo--cannabinoids")
if not cannabinoid_info or "No information provided" in cannabinoid_info.text:
    for x in range(0, 6):
        chemical_content_rows[0].append("n/a")
else:
    # Find all dd elements within cannabinoid section
    cannabinoid_dds = cannabinoid_info.find_all("dd")
    for dd in cannabinoid_dds:
        # Remove % symbol and convert to string
        value = dd.text.strip().replace("%", "")
        chemical_content_rows[0].append(value)

# Extract terpenoid data
terpenoid_info = soup.find(class_="StrainChemicalInfo--terpenoids")
if not terpenoid_info or "No information provided" in terpenoid_info.text:
    for x in range(0, 21):
        chemical_content_rows[0].append("n/a")
else:
    # Find all dd elements within terpenoid section that have values
    terpenoid_dds = terpenoid_info.find_all("dd")
    for dd in terpenoid_dds:
        value = dd.text.strip().replace("%", "")
        if value and value != "n/a":
            chemical_content_rows[0].append(value)
        else:
            chemical_content_rows[0].append("n/a")

# getting variant rows

variants_rows = []

variants = soup.find_all(class_="-js Variants--row")

for item in variants:
    row = []

    # getting fields
    row.append(item.find(attrs={"data-field": "gene"}).text.split()[0])
    row.append(item.find(attrs={"data-field": "hgvsc"}).text.strip())
    row.append(item.find(attrs={"data-field": "hgvsp"}).text.strip())
    row.append(item.find(attrs={"data-field": "annotation"}).text.strip())
    row.append(item.find(attrs={"data-field": "annotation_impact"}).text.strip())
    row.append(item.find(attrs={"data-field": "contig"}).text.strip())
    row.append(item.find(attrs={"data-field": "contig_pos"}).text.split()[0])
    row.append(item.find(attrs={"data-field": "ref_alt"}).text.strip())
    row.append(
        item.find(attrs={"data-field": "var_freq"}).find_all("dd")[0].text.strip()
    )
    row.append(
        item.find(attrs={"data-field": "var_freq"}).find_all("dd")[1].text.strip()
    )

    # adding to the rows of the CSV file
    variants_rows.append(row)

# writing the CSV files: metadata, chemicals, variants

underscored_strain_name = "_".join(strain_name.split())

print("Writing CSV files...")
with open(
    f"{folder_name}/{underscored_strain_name + '.metadata.csv'}", "w", newline=""
) as metadata_csvfile:
    # making a csv writer object
    metadata_csvwriter = csv.writer(metadata_csvfile)

    # writing the fields
    metadata_csvwriter.writerow(metadata_fields)

    # writing the data rows
    metadata_csvwriter.writerows(metadata_rows)

with open(
    f"{folder_name}/{underscored_strain_name + '.chemicals.csv'}", "w", newline=""
) as chemical_csvfile:
    # making a csv writer object
    chemical_csvwriter = csv.writer(chemical_csvfile)

    # writing the fields
    chemical_csvwriter.writerow(chemical_content_fields)

    # writing the data rows
    chemical_csvwriter.writerows(chemical_content_rows)

with open(
    f"{folder_name}/{underscored_strain_name + '.variants.csv'}", "w", newline=""
) as variant_csvfile:
    # making a csv writer object
    variant_csvwriter = csv.writer(variant_csvfile)

    # writing the fields
    variant_csvwriter.writerow(variants_fields)

    # writing the data rows
    variant_csvwriter.writerows(variants_rows)

print("Done!")

def write_readable_summary(folder_name, strain_name, metadata, chemical_data, variants):
    summary_file = os.path.join(folder_name, f"{strain_name}_summary.txt")
    
    with open(summary_file, 'w') as f:
        # Write strain basic info
        f.write(f"{'='*80}\n")
        f.write(f"{strain_name} ({metadata[1]}) Summary\n")
        f.write(f"{'='*80}\n\n")
        
        # General Information
        f.write("GENERAL INFORMATION\n")
        f.write(f"{'-'*80}\n")
        f.write(f"Grower: {metadata[2]}\n")
        f.write(f"Sample Name: {metadata[3]}\n")
        f.write(f"Accession Date: {metadata[4]}\n")
        f.write(f"Reported Sex: {metadata[5]}\n")
        f.write(f"Report Type: {metadata[6]}\n")
        f.write(f"Rarity: {metadata[7]}\n")
        f.write(f"Plant Type: {metadata[8]}\n")
        f.write(f"Reported Heterozygosity: {metadata[9]}\n")
        f.write(f"Y Ratio Distribution: {metadata[10]}\n\n")
        
        # Chemical Content
        f.write("CHEMICAL CONTENT\n")
        f.write(f"{'-'*80}\n")
        f.write("Cannabinoids:\n")
        cannabinoids = [
            ("THC + THCA", chemical_data[0]),
            ("CBD + CBDA", chemical_data[1]),
            ("THCV + THCVA", chemical_data[2]),
            ("CBC + CBCA", chemical_data[3]),
            ("CBG + CBGA", chemical_data[4]),
            ("CBN + CBNA", chemical_data[5])
        ]
        for name, value in cannabinoids:
            if value != "n/a":
                f.write(f"  {name}: {value}%\n")
        
        f.write("\nTerpenoids:\n")
        terpenoids = list(zip(chemical_content_fields[6:], chemical_data[6:]))
        for name, value in terpenoids:
            if value != "n/a":
                f.write(f"  {name}: {value}%\n")
        f.write("\n")
        
        # Genetic Relationships
        f.write("GENETIC RELATIONSHIPS\n")
        f.write(f"{'-'*80}\n")
        
        def format_relationships(relationships_str):
            if not relationships_str:
                return []
            relations = []
            for rel in relationships_str.split("|"):
                if ":" in rel:
                    distance, strain_info = rel.split(":", 1)
                    relations.append((float(distance), strain_info))
            return sorted(relations, key=lambda x: x[0])
        
        # All Samples
        f.write("Nearest Genetic Relatives (All Samples):\n")
        all_samples = format_relationships(metadata[11])
        for distance, strain in all_samples:
            f.write(f"  {distance:.3f} - {strain}\n")
        f.write("\n")
        
        # Base Tree
        f.write("Nearest Genetic Relatives (Base Tree):\n")
        base_tree = format_relationships(metadata[12])
        for distance, strain in base_tree:
            f.write(f"  {distance:.3f} - {strain}\n")
        f.write("\n")
        
        # Most Distant
        f.write("Most Genetically Distant Strains:\n")
        distant = format_relationships(metadata[13])
        for distance, strain in distant:
            f.write(f"  {distance:.3f} - {strain}\n")
        f.write("\n")
        
        # Blockchain Information
        f.write("BLOCKCHAIN INFORMATION\n")
        f.write(f"{'-'*80}\n")
        f.write(f"Transaction ID: {metadata[-2]}\n")
        f.write(f"SHASUM Hash: {metadata[-1]}\n")

write_readable_summary(folder_name, "_".join(strain_name.split()), metadata_rows[0], chemical_content_rows[0], variants_rows)

print(f"Created readable summary in {folder_name}/{strain_name}_summary.txt")