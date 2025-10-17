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
    NCBI_EMAIL = os.getenv('NCBI_EMAIL', 'thomas@edenaglabs.com')
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