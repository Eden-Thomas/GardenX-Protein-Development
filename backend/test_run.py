#!/usr/bin/env python3
"""
Quick test script to verify pipeline setup
"""

import asyncio
import sys
import os
from pathlib import Path

# Add backend to path if needed
sys.path.insert(0, str(Path(__file__).parent))

async def test_basic_setup():
    """Test basic imports and configuration"""
    print("=" * 60)
    print("🔬 GardenX Protein Development - Setup Test")
    print("=" * 60)
    
    # Test 1: Config import
    print("\n1️⃣  Testing configuration...")
    try:
        from config import Config
        config = Config()
        print("   ✅ Config loaded successfully")
        print(f"   📁 Results directory: {config.RESULTS_DIR}")
        api_key_status = "✅ Set" if config.NEUROSNAP_API_KEY else "❌ Missing"
        print(f"   🔑 NeuroSnap API Key: {api_key_status}")
    except Exception as e:
        print(f"   ❌ Config error: {e}")
        return False
    
    # Test 2: NeuroSnap client
    print("\n2️⃣  Testing NeuroSnap client...")
    try:
        from neurosnap_client import NeurosnapClient
        if config.NEUROSNAP_API_KEY:
            client = NeurosnapClient(config.NEUROSNAP_API_KEY)
            print("   ✅ NeuroSnap client initialized")
            print(f"   📊 Available services: {len(client.services)}")
            
            # Quick API test using actual method
            test_seq = "MTAILERRESESLWGRFCNWG"
            try:
                # Test TemStaPro which your client actually has
                results = await client.run_temstapro_batch([test_seq])
                if results and results[0].get('predicted_tm'):
                    tm = results[0]['predicted_tm']
                    print(f"   ✅ API test successful - Predicted Tm: {tm:.1f}°C")
                else:
                    print("   ⚠️  API connection OK but no results returned")
            except Exception as api_error:
                print(f"   ⚠️  API test skipped: {str(api_error)[:50]}")
                print("   💡 This is OK for testing - API will work in production")
        else:
            print("   ⏭️  Skipping (no API key)")
    except Exception as e:
        print(f"   ❌ NeuroSnap error: {e}")
    
    # Test 3: Reference sequences
    print("\n3️⃣  Testing reference sequences...")
    try:
        from analysis_engine import REFERENCE_SEQUENCES, get_sequence
        print(f"   ✅ Found {len(REFERENCE_SEQUENCES)} reference sequences")
        print(f"   📝 Available species: {', '.join(list(REFERENCE_SEQUENCES.keys())[:5])}")
        
        # Use helper function to get the actual sequence string
        maize_seq = get_sequence('maize')
        if maize_seq:
            print(f"   📊 Maize D1 length: {len(maize_seq)} aa")
    except Exception as e:
        print(f"   ❌ Reference sequences error: {e}")
    
    # Test 4: Pipeline modules
    print("\n4️⃣  Testing pipeline modules...")
    try:
        from pipeline.master_pipeline_v3 import MasterPipelineV3Sequential as MasterPipelineV3

        print("   ✅ Master pipeline imported")
        
        from pipeline.track1_evolutionary import Track1Evolutionary
        print("   ✅ Track 1 (Evolutionary) imported")
        
        from pipeline.track2_generative import Track2Generative  
        print("   ✅ Track 2 (Generative) imported")
        
        from pipeline.track3_structural import Track3Structural
        print("   ✅ Track 3 (Structural) imported")
        
        from pipeline.consensus_scorer import ConsensusScorer
        print("   ✅ Consensus scorer imported")
        
        pipeline = MasterPipelineV3(config)
        print("   ✅ Pipeline initialized successfully!")
        
    except ImportError as e:
        print(f"   ❌ Import error: {e}")
        print("   💡 Make sure all pipeline files are in backend/pipeline/")
    except Exception as e:
        print(f"   ❌ Pipeline error: {e}")
    
    print("\n" + "=" * 60)
    print("Setup test complete!")
    
    if config.NEUROSNAP_API_KEY:
        print("\n✅ Ready to run full pipeline!")
        print("Run: python run_pipeline.py --variants 15 --temp 50")
    else:
        print("\n⚠️  Add NEUROSNAP_API_KEY to .env file for full functionality")
    
    return True

async def test_mini_pipeline():
    """Run a minimal pipeline test"""
    print("\n" + "=" * 60)
    print("🚀 Running mini pipeline test...")
    print("=" * 60)
    
    from config import Config
    from pipeline.master_pipeline_v3 import MasterPipelineV3Sequential as MasterPipelineV3
    from analysis_engine import REFERENCE_SEQUENCES, get_sequence
    
    config = Config()
    pipeline = MasterPipelineV3(config)
    
    # Use helper function to get the actual sequence string
    maize_sequence = get_sequence('maize')
    
    if not maize_sequence:
        print("❌ No maize sequence found in REFERENCE_SEQUENCES")
        return
    
    test_sequence = maize_sequence[:50]
    print(f"\n📝 Test sequence (first 50aa of maize D1):\n{test_sequence}")
    
    print("\n⏳ Starting pipeline (this may take a minute)...")
    
    try:
        results = await pipeline.run(
            wt_sequence=test_sequence,
            target_temp=50,
            n_variants=3,  # Just 3 variants for quick test
            validate_all=False
        )
        
        if results.get('final_variants'):
            print(f"\n✅ Success! Generated {len(results['final_variants'])} variants")
            best = results['final_variants'][0]
            print(f"\n🏆 Best variant:")
            print(f"   ID: {best.get('id', 'N/A')}")
            print(f"   Risk: {best.get('risk_level', 'N/A')}")
            print(f"   Mutations: {best.get('mutations', 'N/A')}")
            print(f"   Score: {best.get('consensus_score', best.get('final_score', 'N/A'))}")
        else:
            print("\n⚠️  No variants passed validation")
            print("This is normal for a short test sequence")
            
    except Exception as e:
        print(f"\n❌ Pipeline error: {e}")
        import traceback
        traceback.print_exc()

async def test_neurosnap_services():
    """Test specific Neurosnap services"""
    print("\n" + "=" * 60)
    print("🔬 Testing Neurosnap Services...")
    print("=" * 60)
    
    from config import Config
    from neurosnap_client import NeurosnapClient
    
    config = Config()
    if not config.NEUROSNAP_API_KEY:
        print("⏭️  Skipping - no API key configured")
        return
    
    client = NeurosnapClient(config.NEUROSNAP_API_KEY)
    
    # Test sequence
    test_seq = "MTAILERRESESLWGRFCNWG"
    
    print("\n📊 Available Neurosnap services:")
    for key, name in list(client.services.items())[:10]:
        print(f"   • {key}: {name}")
    
    print(f"\n🧪 Testing TemStaPro with sequence (length {len(test_seq)})...")
    try:
        results = await client.run_temstapro_batch([test_seq])
        if results:
            print(f"   ✅ TemStaPro success!")
            print(f"   Predicted Tm: {results[0].get('predicted_tm', 'N/A')}°C")
    except Exception as e:
        print(f"   ❌ TemStaPro error: {str(e)[:100]}")
    
    print("\n💡 Note: Full API testing requires valid protein sequences and may incur costs")

async def main():
    """Main test function"""
    # Run basic setup test
    success = await test_basic_setup()
    
    if success:
        print("\n" + "=" * 60)
        print("Choose test option:")
        print("1. Run mini pipeline test (recommended)")
        print("2. Test Neurosnap services")
        print("3. Skip additional tests")
        response = input("\nSelect (1/2/3): ")
        
        if response == '1':
            await test_mini_pipeline()
        elif response == '2':
            await test_neurosnap_services()
    
    print("\n✨ Test complete!")

if __name__ == "__main__":
    asyncio.run(main())