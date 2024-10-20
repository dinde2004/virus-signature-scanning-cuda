import re
import sys

# Check if a file path is provided
if len(sys.argv) != 2:
    print("Usage: python script.py <file_path>")
    sys.exit(1)

# Reading the file path from the command-line argument
file_path = sys.argv[1]

# Read the contents of the file
try:
    with open(file_path, 'r') as file:
        input_text = file.read()
except FileNotFoundError:
    print(f"Error: File {file_path} not found.")
    sys.exit(1)

# Regular expression to match the time
time_pattern = r"Total runMatcher time:(\d+\.\d+)s"

# Find all matches in the text
times = re.findall(time_pattern, input_text)

# Print each extracted time on a new line
for time in times:
    print(time)
