import requests

# URL of the file to download
url = "http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2020-02-25.txt"

# Name of the output file to save the downloaded data
output_filename = "DGV_data.csv"

# Send a GET request to the URL
response = requests.get(url)

# Check if the request was successful (status code 200)
if response.status_code == 200:
    # Open the output file in binary write mode
    with open(output_filename, 'wb') as file:
        # Write the content of the response (downloaded data) to the file
        file.write(response.content)
    # Print a message indicating that the file has been downloaded and saved
    print(f"File downloaded and saved as '{output_filename}'")
else:
    # If the request was not successful, print an error message
    print("Failed to download the file")
