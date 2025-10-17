"""
run_pipeline.py
Main execution script for D1 engineering pipeline
"""

import asyncio
import os
from dotenv import load_dotenv
from pipeline.master_pipeline_v3 import D1EngineeringMasterPipeline
from analysis_engine import REFERENCE_SEQUENCES

load_dotenv()

async def main():
    # Get API key
    api_key = os.getenv('NEUROSNAP_API_KEY')
    if not api_key:
        print("❌ Please set NEUROSNAP_API_KEY in .env file")
        return
    
    # Initialize pipeline
    pipeline = D1EngineeringMasterPipeline(api_key)
    
    # Use corn D1 as wild-type
    wt_sequence = REFERENCE_SEQUENCES['corn']
    
    # Run pipeline
    results = await pipeline.run(
        wt_sequence=wt_sequence,
        target_temp=50,  # Target 50°C (vs ~28°C native)
        n_variants=100,  # Start small for testing
        top_k=10  # Top 10 final variants
    )
    
    # Print summary
    print("\n" + "="*50)
    print("FINAL RESULTS")
    print("="*50)
    
    for i, variant in enumerate(results['variants'], 1):
        print(f"\n{i}. {variant['id']}")
        print(f"   Composite Score: {variant['validation']['composite_score']:.1f}/100")
        print(f"   ΔTm: {variant['validation']['delta_tm']:+.1f}°C")
        print(f"   Risk: {variant['risk_level']}")
        print(f"   Mutations: {len(variant.get('mutations', []))}")

if __name__ == "__main__":
    asyncio.run(main())