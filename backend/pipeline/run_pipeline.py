"""
run_pipeline.py
Main execution script for D1 engineering pipeline
Supports both parallel and sequential execution modes
"""

import asyncio
import argparse
import sys
import os
from pathlib import Path
from dotenv import load_dotenv
from datetime import datetime

# Add parent directory to path
sys.path.append(str(Path(__file__).parent))

# Load environment variables
load_dotenv()

# Import configuration and utilities
from config import Config
from analysis_engine import get_sequence

async def main():
    """Main execution function for D1 protein engineering pipeline"""
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='D1 Protein Engineering Pipeline - Generate heat-tolerant variants',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate 50 variants for 50°C using sequential execution
  python run_pipeline.py --variants 50 --temp 50
  
  # Use sequential pipeline where each track builds on previous insights
  python run_pipeline.py --variants 100 --temp 55 --sequential
  
  # Run with full validation and maize sequence
  python run_pipeline.py --sequence maize --validate-all --sequential
        """
    )
    
    parser.add_argument('--variants', type=int, default=50, 
                       help='Number of variants to generate (default: 50)')
    parser.add_argument('--temp', type=int, default=50,
                       help='Target temperature in Celsius (default: 50)')
    parser.add_argument('--sequence', type=str, default='maize',
                       choices=['maize', 'sorghum', 'wheat', 'rice', 'setaria', 'panicum'],
                       help='Reference sequence to use (default: maize)')
    parser.add_argument('--validate-all', action='store_true',
                       help='Run full validation on all variants (slower but more thorough)')
    parser.add_argument('--sequential', action='store_true', default=True,
                       help='Use sequential pipeline where Track1→Track2→Track3 with insights passed between (DEFAULT)')
    parser.add_argument('--output-dir', type=str, default=None,
                       help='Custom output directory for results')
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug output')
    
    args = parser.parse_args()
    
    # Print header
    print("\n" + "="*70)
    print("    D1 PROTEIN ENGINEERING PIPELINE")
    print("    Heat Tolerance Optimization for Photosystem II")
    print("="*70)
    
    # Initialize configuration
    config = Config()
    
    # Override config with command line arguments
    config.TARGET_TEMPERATURE = args.temp
    config.INITIAL_VARIANT_COUNT = args.variants
    
    if args.output_dir:
        config.RESULTS_DIR = Path(args.output_dir)
        config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # FIXED: Always use sequential pipeline (Track 1 → Track 2 → Track 3)
    # This is the correct mode where each track builds on previous insights
    try:
        from pipeline.master_pipeline_v3 import MasterPipelineV3Sequential
        pipeline = MasterPipelineV3Sequential(config)
        print("🔗 Mode: SEQUENTIAL (Track 1 → Track 2 → Track 3)")
        print("   Each track builds on previous insights")
        print("   Full statistical analysis enabled")
    except ImportError as e:
        print(f"❌ Error: Sequential pipeline not found")
        print(f"   Error details: {e}")
        print(f"   Make sure pipeline/master_pipeline_v3.py has MasterPipelineV3Sequential class")
        return
    
    # Get wild-type sequence
    wt_sequence = get_sequence(args.sequence)
    
    if not wt_sequence:
        print(f"\n❌ Error: Could not find sequence for '{args.sequence}'")
        print("Available sequences: maize, sorghum, wheat, rice, setaria, panicum")
        return
    
    # Display run parameters
    print(f"\n📊 Input Parameters:")
    print(f"   • Reference: {args.sequence} D1 protein ({len(wt_sequence)} aa)")
    print(f"   • Target temperature: {args.temp}°C")
    print(f"   • Variants to generate: {args.variants}")
    print(f"   • Full validation: {'Yes' if args.validate_all else 'No'}")
    print(f"   • Output directory: {config.RESULTS_DIR}")
    
    # Check API key
    if not config.NEUROSNAP_API_KEY:
        print("\n⚠️  Warning: No NEUROSNAP_API_KEY found")
        print("   Pipeline will run with limited functionality")
        response = input("   Continue anyway? (y/n): ")
        if response.lower() != 'y':
            print("Exiting...")
            return
    else:
        print(f"\n✅ API Key configured")
    
    print("\n" + "="*70)
    print("Starting pipeline execution...")
    print("="*70 + "\n")
    
    # Record start time
    start_time = datetime.now()
    
    try:
        # Run the pipeline
        results = await pipeline.run(
            wt_sequence=wt_sequence,
            target_temp=args.temp,
            n_variants=args.variants,
            validate_all=args.validate_all
        )
        
        # Calculate runtime
        runtime = (datetime.now() - start_time).total_seconds() / 60
        
        # Process results
        n_final = len(results.get('final_variants', []))
        
        print("\n" + "="*70)
        print("    PIPELINE EXECUTION COMPLETE")
        print("="*70)
        print(f"\n📊 Final Results Summary:")
        print(f"   • Total runtime: {runtime:.1f} minutes")
        print(f"   • Validated variants: {n_final}")
        
        if n_final > 0:
            # Show top 3 variants
            print(f"\n🏆 Top 3 Variants:")
            print("─" * 50)
            
            for i, variant in enumerate(results['final_variants'][:3], 1):
                print(f"\n  {i}. Variant {variant.get('id', 'Unknown')}")
                print(f"     Mutations: {', '.join(variant.get('mutations', []))}")
                print(f"     Predicted ΔTm: +{variant.get('delta_tm', 0):.1f}°C")
                print(f"     Risk level: {variant.get('risk_level', 'Unknown')}")
                print(f"     Confidence: {variant.get('confidence_score', 0):.1%}")
                
                # Show mechanistic rationale if using sequential pipeline
                if args.sequential and 'rationale' in variant:
                    rationale = variant['rationale']
                    mechanisms = rationale.get('mechanisms', [])
                    if mechanisms:
                        print(f"     Why it works:")
                        for mechanism in mechanisms[:2]:
                            print(f"       • {mechanism}")
                    
                    stats = rationale.get('statistical_support', {})
                    if stats:
                        print(f"     Statistical support: {stats.get('p_value', 'N/A')}")
            
            # Show statistical summary if available
            if args.sequential and 'statistics' in results:
                stats = results['statistics'].get('overall_stats', {})
                if stats:
                    print(f"\n📊 Statistical Analysis:")
                    print(f"   • Mean predicted Tm: {stats.get('mean_tm', 0):.1f}°C")
                    print(f"   • Best ΔTm achieved: +{stats.get('best_delta_tm', 0):.1f}°C")
                    print(f"   • Statistical significance: p<{stats.get('significance', 0.05):.4f}")
                    print(f"   • Correlation (mutations vs Tm): r={stats.get('correlation', 0):.3f}")
            
            # Output file locations
            print(f"\n📁 Output Files:")
            output_dir = config.RESULTS_DIR
            print(f"   • Directory: {output_dir}")
            
            # List generated files
            if output_dir.exists():
                files = list(output_dir.glob(f"*{datetime.now().strftime('%Y%m%d')}*"))
                if files:
                    print(f"   • Generated files:")
                    for f in sorted(files)[-5:]:  # Show last 5 files
                        print(f"     - {f.name}")
        
        else:
            print("\n⚠️  No variants passed validation")
            print("   This may be due to:")
            print("   • SSL/connection issues with the API")
            print("   • Very stringent validation criteria")
            print("   • Short test sequence")
            
            if args.debug:
                print("\n🔍 Debug Information:")
                if 'error' in results:
                    print(f"   Error: {results['error']}")
                if 'statistics' in results:
                    print(f"   Pipeline stats: {results['statistics']}")
        
        print("\n" + "="*70)
        print("✅ Pipeline execution finished successfully")
        print("="*70)
        
    except KeyboardInterrupt:
        print("\n\n⚠️  Pipeline interrupted by user")
        print("Partial results may have been saved to:", config.RESULTS_DIR)
        
    except Exception as e:
        print(f"\n❌ Pipeline error: {e}")
        if args.debug:
            import traceback
            print("\n🔍 Full traceback:")
            traceback.print_exc()
        print("\nTroubleshooting tips:")
        print("  1. Check your NEUROSNAP_API_KEY in .env")
        print("  2. Verify SSL/TLS configuration")
        print("  3. Try with fewer variants (--variants 10)")
        print("  4. Enable debug mode (--debug)")


if __name__ == "__main__":
    # Check Python version
    if sys.version_info < (3, 8):
        print("❌ Python 3.8+ is required")
        sys.exit(1)
    
    # Run the async main function
    asyncio.run(main())