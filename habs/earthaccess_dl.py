import earthaccess
import os

# 1. AUTHENTICATION
# This will open a link or prompt for your NASA Earthdata credentials.
# It saves them to your .netrc, so you only do this once.
earthaccess.login(strategy="interactive")

# 2. DEFINE SEARCH PARAMETERS
# We are looking for "VIIRS SNPP Level-3 Mapped Chlorophyll-a" 
# This product is already projected to a 2D lat/lon grid.
DATASET_SHORT_NAME = "VIIRSN_L3m_CHL" 

# Lake Erie Bounding Box: [Lower-Left Lon, Lat, Upper-Right Lon, Lat]
ERIE_BBOX = (-83.6, 41.3, -78.2, 43.1)

# Set your date range (e.g., the August HAB season)
DATE_RANGE = ("2025-08-01", "2025-08-31")

# 3. SEARCH
print(f"Searching for {DATASET_SHORT_NAME} over Lake Erie...")
results = earthaccess.search_data(
    short_name=DATASET_SHORT_NAME,
    bounding_box=ERIE_BBOX,
    temporal=DATE_RANGE
)

print(f"Found {len(results)} files matching your criteria.")

# 4. DOWNLOAD
# Create a local directory for the data
output_folder = "./lake_erie_viirs_data"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

if len(results) > 0:
    print(f"Downloading files to {output_folder}...")
    # This downloads the files in parallel for speed
    downloaded_files = earthaccess.download(results, output_folder)
    print("\nDownload Complete!")
else:
    print("No data found. Try expanding your date range.")