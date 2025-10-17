"""
setup_apis.py
Interactive setup for NeuroSnap and Evo 2 API keys
"""

from pathlib import Path
import os

def setup_api_keys():
    """Interactive API key setup"""
    
    print("="*70)
    print("API KEY SETUP")
    print("="*70)
    
    print("\nYou need two API keys:")
    print("1. NeuroSnap (for AlphaFold2 + NeuroFold)")
    print("   Get it at: https://neurosnap.ai → Sign up → API tab")
    print("\n2. NVIDIA (for Evo 2)")
    print("   Get it at: https://build.nvidia.com → Evo 2 → Get API key")
    
    print("\n" + "-"*70)
    
    # Get API keys
    print("\nEnter your API keys (or press Enter to skip):")
    neurosnap_key = input("NeuroSnap API Key: ").strip()
    nvidia_key = input("NVIDIA API Key: ").strip()
    
    # Create .env file
    env_file = Path(".env")
    
    env_content = f"""# GardenX D1 Platform - API Keys
# DO NOT COMMIT THIS FILE TO GIT

# NeuroSnap API (AlphaFold2, NeuroFold)
NEUROSNAP_API_KEY={neurosnap_key if neurosnap_key else 'your_neurosnap_key_here'}

# NVIDIA API (Evo 2)
NVIDIA_API_KEY={nvidia_key if nvidia_key else 'your_nvidia_key_here'}
EVO2_API_URL=https://health.api.nvidia.com/v1/biology/arc/evo2-40b
"""
    
    with open(env_file, 'w') as f:
        f.write(env_content)
    
    print(f"\n✓ Created {env_file}")
    
    # Update .gitignore
    gitignore = Path(".gitignore")
    
    gitignore_entries = [".env", "*.env", ".env.local"]
    
    if gitignore.exists():
        with open(gitignore, 'r') as f:
            content = f.read()
        
        # Add .env if not present
        if '.env' not in content:
            with open(gitignore, 'a') as f:
                f.write('\n# Environment variables\n')
                for entry in gitignore_entries:
                    f.write(f'{entry}\n')
            print(f"✓ Updated {gitignore}")
    else:
        with open(gitignore, 'w') as f:
            f.write('# Environment variables\n')
            for entry in gitignore_entries:
                f.write(f'{entry}\n')
        print(f"✓ Created {gitignore}")
    
    # Test if keys work
    print("\n" + "-"*70)
    print("TESTING API CONNECTIONS")
    print("-"*70)
    
    if neurosnap_key:
        test_neurosnap(neurosnap_key)
    else:
        print("\n⚠  NeuroSnap key not provided - skipping test")
    
    if nvidia_key:
        test_nvidia(nvidia_key)
    else:
        print("\n⚠  NVIDIA key not provided - skipping test")
    
    print("\n" + "="*70)
    print("API SETUP COMPLETE!")
    print("="*70)
    
    if not neurosnap_key or not nvidia_key:
        print("\n⚠  Note: You can run the platform in local-only mode")
        print("   without API keys, but you won't get structure validation.")
        print("\n   To add keys later, edit the .env file")


def test_neurosnap(api_key):
    """Test NeuroSnap API connection"""
    import requests
    
    try:
        url = "https://neurosnap.ai/api/services"
        headers = {"X-API-KEY": api_key}
        response = requests.get(url, headers=headers, timeout=10)
        
        if response.status_code == 200:
            print("\n✓ NeuroSnap API: Connected successfully")
            services = response.json()
            print(f"  Available services: {len(services)}")
        else:
            print(f"\n✗ NeuroSnap API: Error {response.status_code}")
            print("  Please check your API key")
    except Exception as e:
        print(f"\n✗ NeuroSnap API: Connection failed - {e}")


def test_nvidia(api_key):
    """Test NVIDIA API connection"""
    import requests
    
    try:
        # Simple test - check if key format is valid
        if api_key and len(api_key) > 10:
            print("\n✓ NVIDIA API: Key format looks valid")
            print("  (Full test requires sequence submission)")
        else:
            print("\n✗ NVIDIA API: Key format invalid")
    except Exception as e:
        print(f"\n✗ NVIDIA API: Error - {e}")


def main():
    """Run API setup"""
    setup_api_keys()
    
    print("\nNext steps:")
    print("1. Download datasets: python backend/download_datasets.py")
    print("2. Integrate data: python backend/integrate_sequences.py")
    print("3. Run analysis: python backend/master_pipeline_v1.py")


if __name__ == "__main__":
    main()