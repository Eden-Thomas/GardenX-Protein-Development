#!/usr/bin/env python3
"""
Debug .env file to check for hidden characters or format issues
"""

import os
import binascii

def debug_env_file():
    """Check for common .env file issues"""
    
    print("=" * 60)
    print("DEBUGGING .ENV FILE")
    print("=" * 60)
    
    # 1. Check if .env exists
    if not os.path.exists('.env'):
        print("❌ .env file not found!")
        return
    
    print("✅ .env file exists\n")
    
    # 2. Read the raw file
    with open('.env', 'rb') as f:
        raw_content = f.read()
    
    # 3. Find the NEUROSNAP_API_KEY line
    lines = raw_content.decode('utf-8', errors='replace').split('\n')
    api_key_line = None
    line_num = 0
    
    for i, line in enumerate(lines, 1):
        if 'NEUROSNAP_API_KEY' in line:
            api_key_line = line
            line_num = i
            break
    
    if not api_key_line:
        print("❌ NEUROSNAP_API_KEY not found in .env file!")
        return
    
    print(f"Found NEUROSNAP_API_KEY on line {line_num}")
    print(f"Raw line: {repr(api_key_line)}")
    print()
    
    # 4. Check for common issues
    issues = []
    
    # Check for BOM (Byte Order Mark)
    if raw_content.startswith(b'\xef\xbb\xbf'):
        issues.append("File has UTF-8 BOM (Byte Order Mark)")
    
    # Check for quotes
    if '"' in api_key_line or "'" in api_key_line:
        issues.append("Line contains quotes (should not have any)")
    
    # Check for spaces around =
    if ' =' in api_key_line or '= ' in api_key_line:
        issues.append("Spaces around = sign (should be NEUROSNAP_API_KEY=key)")
    
    # Check for trailing whitespace
    if api_key_line.rstrip() != api_key_line:
        issues.append(f"Trailing whitespace detected: {repr(api_key_line[-5:])}")
    
    # Check for carriage returns (Windows line endings)
    if '\r' in api_key_line:
        issues.append("Windows line endings (\\r\\n) detected")
    
    # Extract the key value
    if '=' in api_key_line:
        key_part = api_key_line.split('=', 1)[1].strip()
        print(f"Extracted key: {key_part[:20]}...{key_part[-20:]}")
        print(f"Key length: {len(key_part)}")
        
        # Check if it's hex
        try:
            binascii.unhexlify(key_part)
            print("✅ Key is valid hexadecimal")
        except:
            issues.append("Key is not valid hexadecimal")
        
        # Check for non-ASCII characters
        if not all(ord(c) < 128 for c in key_part):
            issues.append("Non-ASCII characters detected in key")
    
    # 5. Report issues
    if issues:
        print("\n⚠️  ISSUES FOUND:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print("\n✅ No formatting issues detected")
    
    # 6. Show the correct format
    print("\n" + "=" * 60)
    print("CORRECT FORMAT:")
    print("NEUROSNAP_API_KEY=555fda54ce175776e8ed1ce95c97466406522d0e2f6f55dc1637f5a3f8baaa93f38170023ee517e0a2b196c704a750e5263a30fa763705b")
    print("(No quotes, no spaces, all on one line)")
    print("=" * 60)
    
    # 7. Load with python-dotenv to see what it gets
    from dotenv import load_dotenv
    load_dotenv(override=True)
    loaded_key = os.getenv('NEUROSNAP_API_KEY')
    
    print(f"\nWhat python-dotenv loaded:")
    if loaded_key:
        print(f"  Key: {loaded_key[:20]}...{loaded_key[-20:]}")
        print(f"  Length: {len(loaded_key)}")
    else:
        print("  ❌ Key not loaded by python-dotenv!")

if __name__ == "__main__":
    debug_env_file()