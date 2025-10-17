#!/usr/bin/env python3
"""
Test Neurosnap API authentication with correct endpoint
"""

from dotenv import load_dotenv
import os
import requests
import json

# Force reload environment variables
load_dotenv(override=True)

# Get the API key
api_key = os.getenv('NEUROSNAP_API_KEY')

if not api_key:
    print("❌ No API key found in environment")
    exit(1)

print(f"Testing with API key: {api_key[:20]}...{api_key[-20:]}")
print(f"Key length: {len(api_key)}")

# Test 1: Try to list available services (this is a valid endpoint)
print("\n1. Testing service list endpoint...")
response = requests.get(
    "https://neurosnap.ai/api/services",
    headers={"X-API-KEY": api_key}
)

print(f"Status: {response.status_code}")
if response.status_code == 200:
    services = response.json()
    print(f"✅ Authentication successful! Found {len(services)} services")
    # Show first few services
    for service in services[:3]:
        print(f"  - {service}")
else:
    print(f"❌ Error: {response.text}")

# Test 2: Try to list jobs (another valid endpoint)
print("\n2. Testing jobs list endpoint...")
response = requests.get(
    "https://neurosnap.ai/api/jobs",
    headers={"X-API-KEY": api_key}
)

print(f"Status: {response.status_code}")
if response.status_code == 200:
    jobs = response.json()
    print(f"✅ Can list jobs! Found {len(jobs)} jobs")
else:
    print(f"❌ Error: {response.text}")

# Test 3: Test the exact format your neurosnap_client uses
print("\n3. Testing with session (like your client)...")
session = requests.Session()
session.headers.update({"X-API-KEY": api_key})

response = session.get("https://neurosnap.ai/api/services")
print(f"Status: {response.status_code}")
if response.status_code == 200:
    print("✅ Session-based auth works!")
else:
    print(f"❌ Session error: {response.text}")

# Test 4: Verify the key format is correct
print("\n4. Checking API key format...")
if len(api_key) == 128 and api_key.isalnum():
    print("✅ API key format looks correct (128 hex characters)")
else:
    print(f"⚠️  API key format may be wrong. Length: {len(api_key)}, Alphanumeric: {api_key.isalnum()}")