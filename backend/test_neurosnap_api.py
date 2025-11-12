#!/usr/bin/env python3
"""
Test script to verify Neurosnap API connectivity and basic functionality
Eden Agriculture - Testing real API endpoints
"""

import asyncio
import os
import json
from pathlib import Path

# Try to load .env file if it exists
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

from neurosnap_client import NeurosnapAPIClient, D1_SEQUENCES

async def test_api_connection():
    """Test basic API connectivity and endpoints"""
    
    # Set debug mode
    os.environ["DEBUG"] = "true"
    
    # API Key
    API_KEY = os.getenv("NEUROSNAP_API_KEY", "2ac33d9b432efc48bd2ba2aeaaaf260d990f7fe1b2f85d8fc36ab694574fde7975f056c5d76a6399717443ab09355984fcb607aa345e41b80c700ea710eb8c66")
    
    print("\n" + "="*60)
    print("NEUROSNAP API CONNECTION TEST")
    print("="*60 + "\n")
    
    async with NeurosnapAPIClient(API_KEY) as client:
        
        # Test 1: Fitness prediction
        print("Test 1: Fitness Prediction")
        print("-" * 30)
        test_sequence = D1_SEQUENCES["wheat"][:100]  # Use first 100 amino acids for testing
        
        try:
            result = await client.predict_fitness([test_sequence])
            if "error" in result:
                print(f"❌ Fitness prediction failed: {result['error']}")
            else:
                print(f"✅ Fitness prediction successful")
                print(f"   Response keys: {list(result.keys())}")
                if result.get("predictions"):
                    print(f"   Predictions: {json.dumps(result['predictions'], indent=2)[:200]}...")
        except Exception as e:
            print(f"❌ Error: {e}")
        
        print()
        
        # Test 2: Structure prediction
        print("Test 2: Structure Prediction")
        print("-" * 30)
        
        try:
            result = await client.predict_structure(test_sequence[:50])  # Even shorter for structure
            if "error" in result:
                print(f"❌ Structure prediction failed: {result.get('error', 'Unknown error')}")
            else:
                print(f"✅ Structure prediction successful")
                print(f"   Mean pLDDT: {result.get('mean_plddt', 0):.2f}")
                print(f"   PDB available: {bool(result.get('pdb_string'))}")
        except Exception as e:
            print(f"❌ Error: {e}")
        
        print()
        
        # Test 3: Variant generation
        print("Test 3: Variant Generation")
        print("-" * 30)
        
        try:
            # Generate variants for a small region
            target_positions = [105, 110, 115]  # Test positions
            result = await client.generate_variants(
                template_sequence=D1_SEQUENCES["wheat"],
                target_positions=target_positions,
                n_variants=2,  # Just 2 variants for testing
                method="progen2"
            )
            
            if isinstance(result, dict) and "error" in result:
                print(f"❌ Variant generation failed: {result['error']}")
            else:
                print(f"✅ Variant generation successful")
                print(f"   Generated {len(result)} variants")
                if result:
                    print(f"   First variant mutations: {result[0].get('mutations', [])}")
        except Exception as e:
            print(f"❌ Error: {e}")
        
        print()
        
        # Test 4: Thermostability prediction
        print("Test 4: Thermostability Prediction")
        print("-" * 30)
        
        try:
            result = await client.compute_thermostability([test_sequence])
            if "error" in result:
                print(f"❌ Thermostability prediction failed: {result['error']}")
            else:
                print(f"✅ Thermostability prediction successful")
                print(f"   Response keys: {list(result.keys())}")
        except Exception as e:
            print(f"❌ Error: {e}")
    
    print("\n" + "="*60)
    print("API CONNECTIVITY TEST COMPLETE")
    print("="*60 + "\n")

if __name__ == "__main__":
    asyncio.run(test_api_connection())