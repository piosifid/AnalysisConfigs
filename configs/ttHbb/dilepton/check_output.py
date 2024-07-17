from coffea.util import load
import sys
import collections

if len(sys.argv) != 2:
    print("Usage: python check_empty_keys.py <path_to_your_file>")
    sys.exit(1)

filename = sys.argv[1]

try:
    data = load(filename)
    
    # List of dictionaries to inspect
    dictionaries_to_inspect = [
        'variables',
        'columns',
        'processing_metadata',
        'datasets_metadata',
        'trigger_efficiency'  # Add more dictionaries as needed
    ]
    
    # Iterate through each dictionary and check if keys are empty
    for dictionary_name in dictionaries_to_inspect:
        if dictionary_name in data:
            dictionary_data = data[dictionary_name]
            
            print(f"Checking '{dictionary_name}' in '{filename}':")
            
            if isinstance(dictionary_data, dict) and len(dictionary_data) > 0:
                print(f"  '{dictionary_name}' is not empty.")
                # Optionally, print number of keys and some summary of each key if needed
                keys = list(dictionary_data.keys())
                print(f"  Keys found: {keys}")
                # Example: Print number of entries in each key
                for key in keys:
                    print(f"    Number of entries in '{key}': {len(dictionary_data[key])}")
            else:
                print(f"  '{dictionary_name}' is empty.")
            
            print()  # Add a blank line for separation
        else:
            print(f"No '{dictionary_name}' found in '{filename}'.")
            print()  # Add a blank line for separation

except Exception as e:
    print(f"Error loading {filename}: {e}")

