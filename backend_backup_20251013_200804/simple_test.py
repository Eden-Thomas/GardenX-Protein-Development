import sys
print("Python version:", sys.version)
print("Starting import test...")

try:
    from neurosnap_client import NeurosnapClient
    print("✅ NeurosnapClient imported")
except Exception as e:
    print(f"❌ Import error: {e}")

try:
    import aiohttp
    print("✅ aiohttp imported")
except Exception as e:
    print(f"❌ aiohttp not installed: {e}")

try:
    from dotenv import load_dotenv
    print("✅ dotenv imported")
    load_dotenv()
except Exception as e:
    print(f"❌ dotenv error: {e}")

import os
api_key = os.getenv('NEUROSNAP_API_KEY')
if api_key:
    print(f"✅ API key found: {api_key[:20]}...")
else:
    print("❌ No API key in environment")

print("\nAll checks complete!")
