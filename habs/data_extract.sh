#!/bin/bash

# Define directories
INPUT_DIR="./zips"
OUTPUT_DIR="./extracted_chl"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all zip files in the input directory
for zip_file in "$INPUT_DIR"/*.zip; do
    # Check if files exist to avoid errors on empty directories
    [ -e "$zip_file" ] || continue

    # 1. Extract the timestamp from the filename using grep
    # This looks for the 8digit+T+6digit pattern (e.g., 20230815T153826)
    TIMESTAMP=$(echo "$zip_file" | grep -oE '[0-9]{8}T[0-9]{6}' | head -1)

    if [ -z "$TIMESTAMP" ]; then
        echo "Skipping $zip_file: No timestamp found."
        continue
    fi

    # 2. Identify the internal path of chl_nn.nc
    # Sentinel-3 zips contain a .SEN3 folder, so we need the full internal path
    INTERNAL_PATH=$(unzip -l "$zip_file" | grep -oE '[^ ]*chl_nn\.nc' | head -1)

    if [ -n "$INTERNAL_PATH" ]; then
        # 3. Extract and rename
        # -p extracts to pipe (stdout), which we redirect to our new filename
        unzip -p "$zip_file" "$INTERNAL_PATH" > "$OUTPUT_DIR/chl_nn_${TIMESTAMP}.nc"
        echo "Extracted: chl_nn_${TIMESTAMP}.nc"
    else
        echo "Error: chl_nn.nc not found in $zip_file"
    fi

    # 4. Identify the internal path of geo_coordinates.nc
    INTERNAL_PATH=$(unzip -l "$zip_file" | grep -oE '[^ ]*geo_coordinates\.nc' | head -1)
    if [ -n "$INTERNAL_PATH" ]; then
        unzip -p "$zip_file" "$INTERNAL_PATH" > "$OUTPUT_DIR/geo_coordinates_${TIMESTAMP}.nc"
        echo "Extracted: geo_coordinates_${TIMESTAMP}.nc"
    else
        echo "Error: geo_coordinates.nc not found in $zip_file"
    fi


done

echo "Done!"