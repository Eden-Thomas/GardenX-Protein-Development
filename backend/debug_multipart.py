import requests
import os
from dotenv import load_dotenv
import json

load_dotenv()

api_key = os.getenv('NEUROSNAP_API_KEY')
test_seq = "MTAILERRESESLWGRFCNWG"

print("Testing different request formats...")
print("="*60)

# Test 1: Regular form data (no MultipartEncoder)
print("\nTest 1: Regular form data")
response = requests.post(
    "https://neurosnap.ai/api/job/submit/TemStaPro Protein Thermostability Prediction",
    headers={"X-API-KEY": api_key},
    data={
        "Input Sequences": json.dumps([test_seq]),
        "Temperature Thresholds": json.dumps([40, 45, 50, 55, 60, 65, 70])
    },
    timeout=10
)
print(f"Status: {response.status_code}")
if response.status_code == 200:
    print(f"✅ SUCCESS! Job ID: {response.json().get('job_id')}")
else:
    print(f"Error: {response.text[:200]}")

# Test 2: JSON body instead of form data
print("\nTest 2: JSON body")
response = requests.post(
    "https://neurosnap.ai/api/job/submit/TemStaPro Protein Thermostability Prediction",
    headers={
        "X-API-KEY": api_key,
        "Content-Type": "application/json"
    },
    json={
        "Input Sequences": [test_seq],
        "Temperature Thresholds": [40, 45, 50, 55, 60, 65, 70]
    },
    timeout=10
)
print(f"Status: {response.status_code}")
if response.status_code == 200:
    print(f"✅ SUCCESS! Job ID: {response.json().get('job_id')}")
else:
    print(f"Error: {response.text[:200]}")

# Test 3: Files format (since it's FFSequenceList, might expect file upload)
print("\nTest 3: Files format")
response = requests.post(
    "https://neurosnap.ai/api/job/submit/TemStaPro Protein Thermostability Prediction",
    headers={"X-API-KEY": api_key},
    files={
        "Input Sequences": ("sequences.json", json.dumps([test_seq]), "application/json"),
        "Temperature Thresholds": ("temps.json", json.dumps([40, 45, 50, 55, 60, 65, 70]), "application/json")
    },
    timeout=10
)
print(f"Status: {response.status_code}")
if response.status_code == 200:
    print(f"✅ SUCCESS! Job ID: {response.json().get('job_id')}")
else:
    print(f"Error: {response.text[:200]}")