import requests
import os
from dotenv import load_dotenv
from requests_toolbelt.multipart.encoder import MultipartEncoder
import json

load_dotenv()

api_key = os.getenv('NEUROSNAP_API_KEY')
test_seq = "MTAILERRESESLWGRFCNWG"

print("Testing different formats for FFSequenceList field type...")
print("=" * 60)

# Based on FFSequenceList, try these formats:
test_formats = [
    # Format 1: List of sequence objects with name and data
    [{"name": "sequence_1", "data": test_seq}],
    
    # Format 2: List of dictionaries with sequence key
    [{"sequence": test_seq}],
    
    # Format 3: Simple list of sequences
    [test_seq],
    
    # Format 4: List with FASTA format
    [{"name": "seq1", "data": f">seq1\n{test_seq}"}],
    
    # Format 5: Object with sequences array
    {"sequences": [test_seq]},
    
    # Format 6: Just the sequence as a list (most likely based on 'aa' type)
    [test_seq]
]

for i, test_input in enumerate(test_formats, 1):
    print(f"\nTest {i}: {str(test_input)[:60]}...")
    
    fields = {
        "Input Sequences": json.dumps(test_input),
        "Temperature Thresholds": json.dumps([40, 50, 60, 70])
    }
    
    multipart = MultipartEncoder(fields=fields)
    
    response = requests.post(
        "https://neurosnap.ai/api/job/submit/TemStaPro Protein Thermostability Prediction",
        headers={
            "X-API-KEY": api_key,
            "Content-Type": multipart.content_type
        },
        data=multipart,
        timeout=30
    )
    
    if response.status_code == 200:
        print(f"  ✅ SUCCESS! Format {i} works!")
        print(f"  Job ID: {response.json().get('job_id')}")
        print(f"  Working format: {test_input}")
        break
    else:
        print(f"  ❌ Failed: {response.text[:100]}")