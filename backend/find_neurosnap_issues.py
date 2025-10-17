#!/usr/bin/env python3
"""
Script to identify and report Neurosnap API formatting issues in your backend code.
Run this in your backend directory to find what needs to be fixed.
"""

import os
import re
from pathlib import Path

def analyze_neurosnap_calls(file_path):
    """Analyze a Python file for Neurosnap API formatting issues."""
    issues = []
    
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
        lines = content.split('\n')
    
    # Check if file contains Neurosnap API calls
    if 'neurosnap.ai' not in content.lower() and 'neurosnap' not in content.lower():
        return None
    
    print(f"\n{'='*80}")
    print(f"ANALYZING: {file_path}")
    print(f"{'='*80}")
    
    # Pattern 1: Check for raw string assignments to input fields
    # Looking for patterns like: "Input Sequence": some_variable (without json.dumps)
    input_field_pattern = r'["\']Input\s+\w+["\']\s*:\s*(?!json\.dumps)'
    for i, line in enumerate(lines, 1):
        if re.search(input_field_pattern, line, re.IGNORECASE):
            if 'json.dumps' not in line:
                issues.append({
                    'line': i,
                    'code': line.strip(),
                    'issue': 'Input field may not be using json.dumps()',
                    'fix': 'Wrap value in json.dumps([{"type": "...", "data": value}])'
                })
    
    # Pattern 2: Check for numeric values without quotes in field assignments
    # Looking for: "Number Samples": 100 or "Temperature": 0.5
    numeric_field_pattern = r'["\'](?:Number|Temperature|Seed|Samples?)["\']\s*:\s*(\d+\.?\d*)\s*[,}]'
    for i, line in enumerate(lines, 1):
        matches = re.finditer(numeric_field_pattern, line, re.IGNORECASE)
        for match in matches:
            if not (line[:match.start()].count('"') % 2 == 0):  # Skip if inside a string
                continue
            issues.append({
                'line': i,
                'code': line.strip(),
                'issue': f'Numeric value {match.group(1)} should be a string',
                'fix': f'Change {match.group(1)} to "{match.group(1)}"'
            })
    
    # Pattern 3: Check for incorrect Content-Type header
    content_type_pattern = r'["\']Content-Type["\']\s*:\s*["\']multipart/form-data["\']'
    for i, line in enumerate(lines, 1):
        if re.search(content_type_pattern, line):
            if 'multipart_data.content_type' not in line and 'encoder.content_type' not in line:
                issues.append({
                    'line': i,
                    'code': line.strip(),
                    'issue': 'Content-Type should use multipart_data.content_type',
                    'fix': 'Change to: "Content-Type": multipart_data.content_type'
                })
    
    # Pattern 4: Check if using requests without MultipartEncoder
    if 'requests.post' in content and 'neurosnap' in content.lower():
        if 'MultipartEncoder' not in content and 'FormData' not in content:
            issues.append({
                'line': 0,
                'code': 'File-level issue',
                'issue': 'Not using MultipartEncoder for multipart/form-data',
                'fix': 'Import and use: from requests_toolbelt.multipart.encoder import MultipartEncoder'
            })
    
    # Pattern 5: Check for direct list/dict assignments (not json.dumps)
    list_dict_pattern = r'["\']Input\s+\w+["\']\s*:\s*[\[{]'
    for i, line in enumerate(lines, 1):
        if re.search(list_dict_pattern, line, re.IGNORECASE):
            if 'json.dumps' not in lines[max(0, i-2):min(len(lines), i+2)]:
                issues.append({
                    'line': i,
                    'code': line.strip(),
                    'issue': 'Input field has list/dict literal (needs json.dumps)',
                    'fix': 'Wrap the entire list/dict in json.dumps(...)'
                })
    
    return issues


def generate_fix_template(file_path, issues):
    """Generate a template showing how to fix the issues."""
    print(f"\n{'='*80}")
    print(f"FIXES NEEDED FOR: {file_path}")
    print(f"{'='*80}\n")
    
    for idx, issue in enumerate(issues, 1):
        print(f"Issue #{idx}:")
        print(f"  Line {issue['line']}: {issue['issue']}")
        print(f"  Current code: {issue['code']}")
        print(f"  Fix: {issue['fix']}")
        print()
    
    print("\nCORRECT FORMAT EXAMPLE:")
    print("-" * 80)
    print("""
import json
from requests_toolbelt.multipart.encoder import MultipartEncoder

# Correct way to format sequence data
sequence_data = ">protein1\\nMKLLILAVVFTGSALA..."

# Create fields dictionary with proper formatting
fields = {
    # Sequence inputs - MUST use json.dumps with list of objects
    "Input Sequence": json.dumps([
        {"type": "fasta", "data": sequence_data}
    ]),
    
    # Numeric fields - MUST be strings, not numbers
    "Number Samples": "100",  # NOT 100
    "Temperature": "0.1",     # NOT 0.1
}

# Use MultipartEncoder
multipart_data = MultipartEncoder(fields=fields)

# Make request with correct Content-Type
response = requests.post(
    f"https://neurosnap.ai/api/job/submit/{service_name}",
    headers={
        "X-API-KEY": api_key,
        "Content-Type": multipart_data.content_type,  # NOT just "multipart/form-data"
    },
    data=multipart_data,
)
    """)


def main():
    """Main function to scan backend directory."""
    backend_dir = Path('.')  # Current directory
    
    print("🔍 Scanning for Neurosnap API formatting issues...")
    print(f"📁 Directory: {backend_dir.absolute()}\n")
    
    all_issues = {}
    python_files = list(backend_dir.glob('**/*.py'))
    
    if not python_files:
        print("❌ No Python files found in current directory.")
        print("💡 Make sure you run this script in your backend directory.")
        return
    
    for py_file in python_files:
        try:
            issues = analyze_neurosnap_calls(py_file)
            if issues:
                all_issues[py_file] = issues
        except Exception as e:
            print(f"⚠️  Error analyzing {py_file}: {e}")
    
    if not all_issues:
        print("✅ No Neurosnap API calls found or no obvious formatting issues detected.")
        print("💡 If you're still getting errors, please share your Neurosnap API code.")
        return
    
    # Report findings
    print(f"\n\n{'='*80}")
    print(f"SUMMARY: Found {len(all_issues)} file(s) with potential issues")
    print(f"{'='*80}")
    
    for file_path, issues in all_issues.items():
        generate_fix_template(file_path, issues)
    
    # Generate requirements
    print(f"\n{'='*80}")
    print("REQUIRED DEPENDENCIES")
    print(f"{'='*80}")
    print("""
Make sure you have these installed:
    pip install requests requests-toolbelt

Import statements needed:
    import json
    from requests_toolbelt.multipart.encoder import MultipartEncoder
    """)


if __name__ == "__main__":
    main()