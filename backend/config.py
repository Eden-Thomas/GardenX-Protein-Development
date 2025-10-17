"""
backend/config.py
Production configuration for D1 Engineering Pipeline
"""

import os
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

class Config:
    """Central configuration for the pipeline"""
    
    # API Keys
    NEUROSNAP_API_KEY = os.getenv('NEUROSNAP_API_KEY')
    NCBI_EMAIL = os.getenv('NCBI_EMAIL', 'your_email@domain.com')
    NCBI_API_KEY = os.getenv('NCBI_API_KEY')  # Optional but recommended
    
    # Pipeline Configuration
    TARGET_TEMPERATURE = int(os.getenv('TARGET_TEMPERATURE', 50))  # Target Tm in Celsius
    INITIAL_VARIANT_COUNT = int(os.getenv('INITIAL_VARIANT_COUNT', 15))  # Start small
    PRODUCTION_VARIANT_COUNT = int(os.getenv('PRODUCTION_VARIANT_COUNT', 1000))
    
    # Validation Thresholds
    MIN_TM_IMPROVEMENT = float(os.getenv('MIN_TM_IMPROVEMENT', 10.0))  # Minimum ΔTm
    MIN_PLDDT = float(os.getenv('MIN_PLDDT', 70.0))  # Minimum folding confidence
    MIN_TM_SCORE = float(os.getenv('MIN_TM_SCORE', 0.85))  # Structural similarity
    MAX_AGGREGATION_SCORE = float(os.getenv('MAX_AGGREGATION_SCORE', 0.5))
    MIN_SOLUBILITY = float(os.getenv('MIN_SOLUBILITY', 0.4))
    
    # API Rate Limiting
    MAX_API_CALLS_PER_RUN = int(os.getenv('MAX_API_CALLS_PER_RUN', 100))
    MAX_CONCURRENT_JOBS = int(os.getenv('MAX_CONCURRENT_JOBS', 5))
    
    # Paths
    BASE_DIR = Path(__file__).parent
    RESULTS_DIR = BASE_DIR / 'results'
    DATA_DIR = BASE_DIR / 'data'
    CACHE_DIR = BASE_DIR / 'cache'
    
    # Create directories
    RESULTS_DIR.mkdir(exist_ok=True)
    DATA_DIR.mkdir(exist_ok=True)
    CACHE_DIR.mkdir(exist_ok=True)
    
    # D1-specific Configuration
    D1_SEQUENCE_LENGTH = 353  # Standard D1 protein length
    
    # Active site residues (1-indexed)
    ACTIVE_SITE_RESIDUES = {
        161: 'H',  # His161 - Mn ligand
        170: 'D',  # Asp170 - critical for water splitting
        189: 'E',  # Glu189 - proton transfer
        215: 'R',  # Arg215 - stabilizes cluster
        254: 'Y',  # Tyr254 - electron transfer
        255: 'F',  # Phe255 - hydrophobic packing
        264: 'H',  # His264 - Mn ligand
        271: 'Y',  # Tyr271 - electron transfer
        332: 'H',  # His332 - Mn ligand
        333: 'E',  # Glu333 - proton transfer
        342: 'D',  # Asp342 - Ca ligand
        344: 'A'   # Ala344 - structural
    }
    
    # Transmembrane helices (1-indexed)
    TM_HELICES = [
        (20, 40),   # TM1
        (70, 90),   # TM2
        (160, 180), # TM3
        (220, 240), # TM4
        (290, 310)  # TM5
    ]
    
    # Cost Management
    ESTIMATED_COST_PER_CALL = 0.10  # Estimated $ per API call
    MAX_BUDGET = float(os.getenv('MAX_BUDGET', 100.0))  # Maximum budget in $
    
    def validate(self):
        """Validate configuration"""
        if not self.NEUROSNAP_API_KEY:
            raise ValueError("NEUROSNAP_API_KEY is required in .env file")
        
        if not self.NCBI_EMAIL:
            print("Warning: NCBI_EMAIL not set, using default")
        
        return True


# ============================================================================
# backend/run_pipeline.py
# Main execution script for production pipeline
# ============================================================================

"""
run_pipeline.py
Main execution script for D1 Engineering Pipeline
Usage: python run_pipeline.py [--variants N] [--temp T] [--validate-all]
"""

import asyncio
import argparse
import sys
import os
from datetime import datetime

# Add backend to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from config import Config
from pipeline.master_pipeline_v3 import MasterPipelineV3
from analysis_engine import REFERENCE_SEQUENCES

async def progress_callback(progress: float, message: str):
    """Progress callback for pipeline"""
    bar_length = 50
    filled = int(bar_length * progress)
    bar = '█' * filled + '░' * (bar_length - filled)
    print(f"\r[{bar}] {progress*100:.1f}% - {message}", end='', flush=True)
    if progress >= 1.0:
        print()  # New line at completion

async def main():
    """Main execution function"""
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run D1 Engineering Pipeline')
    parser.add_argument('--variants', type=int, default=15,
                       help='Number of variants to generate (default: 15)')
    parser.add_argument('--temp', type=int, default=50,
                       help='Target temperature in Celsius (default: 50)')
    parser.add_argument('--validate-all', action='store_true',
                       help='Validate all variants (expensive)')
    parser.add_argument('--sequence', type=str, default='corn',
                       help='Reference sequence: corn, wheat, rice, thermo (default: corn)')
    parser.add_argument('--production', action='store_true',
                       help='Run in production mode with more variants')
    
    args = parser.parse_args()
    
    # Load configuration
    config = Config()
    
    try:
        config.validate()
    except ValueError as e:
        print(f"❌ Configuration error: {e}")
        print("Please set up your .env file with required API keys")
        return 1
    
    # Get reference sequence
    if args.sequence in REFERENCE_SEQUENCES:
        wt_sequence = REFERENCE_SEQUENCES[args.sequence]
        print(f"Using {args.sequence} D1 sequence (length: {len(wt_sequence)})")
    else:
        print(f"❌ Unknown sequence: {args.sequence}")
        print(f"Available: {', '.join(REFERENCE_SEQUENCES.keys())}")
        return 1
    
    # Adjust variant count for production
    if args.production:
        n_variants = config.PRODUCTION_VARIANT_COUNT
        print("🚀 Running in PRODUCTION mode")
    else:
        n_variants = args.variants
        print("🔬 Running in VALIDATION mode")
    
    # Initialize pipeline
    pipeline = MasterPipelineV3(config)
    
    # Run pipeline
    try:
        print("\nStarting D1 Engineering Pipeline...")
        print("=" * 60)
        
        results = await pipeline.run(
            wt_sequence=wt_sequence,
            target_temp=args.temp,
            n_variants=n_variants,
            validate_all=args.validate_all,
            progress_callback=progress_callback if not args.validate_all else None
        )
        
        # Check results
        if results.get('status') == 'failed':
            print(f"\n❌ Pipeline failed: {results.get('error')}")
            return 1
        
        # Success
        final_variants = results.get('final_variants', [])
        
        if final_variants:
            print(f"\n✅ SUCCESS! Generated {len(final_variants)} validated variants")
            print(f"📁 Results saved to: {config.RESULTS_DIR}")
            
            # Print best variant
            best = final_variants[0]
            print(f"\n🏆 Best Variant: {best.get('id')}")
            print(f"   Consensus Score: {best.get('consensus_score', 0):.2f}")
            print(f"   ΔTm: {best.get('validation', {}).get('delta_tm', 0):+.1f}°C")
            print(f"   Risk Level: {best.get('risk_level')}")
            
            return 0
        else:
            print("\n⚠️ No variants passed validation")
            print("Consider adjusting thresholds or increasing variant count")
            return 1
            
    except KeyboardInterrupt:
        print("\n\n⚠️ Pipeline interrupted by user")
        return 1
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    # Run async main
    exit_code = asyncio.run(main())
    sys.exit(exit_code)


# ============================================================================
# backend/.env.example
# Example environment file - copy to .env and fill in values
# ============================================================================

"""
# Neurosnap Configuration
NEUROSNAP_API_KEY=your_neurosnap_api_key_here

# NCBI Configuration (for Track 1)
NCBI_EMAIL=your_email@domain.com
NCBI_API_KEY=optional_but_recommended_for_higher_rate_limits

# Pipeline Configuration
TARGET_TEMPERATURE=50
INITIAL_VARIANT_COUNT=15
PRODUCTION_VARIANT_COUNT=1000

# Validation Thresholds
MIN_TM_IMPROVEMENT=10.0
MIN_PLDDT=70.0
MIN_TM_SCORE=0.85
MAX_AGGREGATION_SCORE=0.5
MIN_SOLUBILITY=0.4

# API Management
MAX_API_CALLS_PER_RUN=100
MAX_CONCURRENT_JOBS=5
MAX_BUDGET=100.0
"""