import requests
import json
from requests_toolbelt.multipart.encoder import MultipartEncoder

api_key = "your_key_here"
url = "https://neurosnap.ai/api/job/submit/temstapro"

# Try the simplest possible format
multipart_data = MultipartEncoder(
    fields={
        'sequence': 'MTAILERRESESLWG',  # Very short test
        'temperature': '50'
    }
)

response = requests.post(
    url,
    headers={
        'X-API-KEY': api_key,
        'Content-Type': multipart_data.content_type
    },
    data=multipart_data
)

print(f"Status: {response.status_code}")
print(f"Response: {response.text}")
