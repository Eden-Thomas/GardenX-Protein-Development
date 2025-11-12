"""
Configuration file for Eden Agriculture D1 Protein Engineering Pipeline
"""

import os
from pathlib import Path
from typing import Dict, Any

# Try to load .env file if it exists
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    # dotenv not installed, will use environment variables directly
    pass

# API Configuration
NEUROSNAP_API_CONFIG = {
    "api_key": os.getenv("NEUROSNAP_API_KEY", "2ac33d9b432efc48bd2ba2aeaaaf260d990f7fe1b2f85d8fc36ab694574fde7975f056c5d76a6399717443ab09355984fcb607aa345e41b80c700ea710eb8c66"),
    "base_url": "https://api.neurosnap.ai",
    "fallback_mode": True,  # Set to False when API domain is allowed
    "rate_limit_delay": 0.5,  # seconds between requests
    "timeout": 300,  # request timeout in seconds
}

# Pipeline Configuration
PIPELINE_CONFIG = {
    "target_crops": ["wheat", "corn", "rice"],
    "variants_per_window": 10,  # Number of variants to generate per design window
    "max_variants_to_screen": 100,  # Maximum variants to screen in one batch
    "conservation_threshold": 0.9,  # Minimum conservation ratio for cross-species transfer
}

# D1 Protein Engineering Parameters
ENGINEERING_CONFIG = {
    "min_plddt_score": 70,  # Minimum pLDDT score for structure confidence
    "min_fitness_score": 0.6,  # Minimum fitness score to consider variant
    "target_tm_increase": 10,  # Target melting temperature increase in Celsius
    "max_mutations_per_variant": 5,  # Maximum mutations per variant
    "preserve_functional_sites": True,  # Always preserve QB binding and other critical sites
}

# Thermostability Optimization Parameters
THERMO_CONFIG = {
    "target_operating_temp": 110,  # Target operating temperature in Fahrenheit
    "min_tm_celsius": 55,  # Minimum melting temperature in Celsius
    "hydrophobicity_range": [0.35, 0.55],  # Optimal hydrophobicity range
    "charge_optimization": "slightly_positive",  # Charge distribution strategy
    "avoid_aggregation_prone": True,  # Filter out aggregation-prone sequences
}

# Output Configuration
OUTPUT_CONFIG = {
    "results_dir": Path("./results"),
    "save_intermediate": True,  # Save intermediate results during pipeline
    "export_formats": ["fasta", "json", "csv"],  # Output file formats
    "max_variants_export": 100,  # Maximum variants to export
}

# Model Selection
MODEL_CONFIG = {
    "fitness_model": "esm2_t33_650M",  # ESM model for fitness prediction
    "structure_model": "esmfold",  # Structure prediction model
    "generation_model": "progen2",  # Variant generation model
    "thermostability_model": "thermonet",  # Thermostability prediction model
}

# Logging Configuration
LOGGING_CONFIG = {
    "level": "INFO",  # Set to DEBUG for detailed logging
    "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    "file": "d1_engineering.log",  # Log file path
}

# Validation Thresholds
VALIDATION_CONFIG = {
    "min_overall_score": 0.7,  # Minimum overall score for elite variants
    "qb_binding_required": True,  # Require intact QB binding site
    "electron_transport_min": 0.6,  # Minimum electron transport score
    "functional_activity_min": 70,  # Minimum predicted activity percentage
}

def get_config():
    """Return complete configuration dictionary"""
    return {
        "api": NEUROSNAP_API_CONFIG,
        "pipeline": PIPELINE_CONFIG,
        "engineering": ENGINEERING_CONFIG,
        "thermo": THERMO_CONFIG,
        "output": OUTPUT_CONFIG,
        "models": MODEL_CONFIG,
        "logging": LOGGING_CONFIG,
        "validation": VALIDATION_CONFIG
    }

def update_for_production():
    """Update configuration for production use with real API"""
    NEUROSNAP_API_CONFIG["fallback_mode"] = False
    PIPELINE_CONFIG["variants_per_window"] = 50
    PIPELINE_CONFIG["max_variants_to_screen"] = 500
    print("Configuration updated for production mode with real API access")

def update_for_testing():
    """Update configuration for testing with reduced parameters"""
    NEUROSNAP_API_CONFIG["fallback_mode"] = True
    PIPELINE_CONFIG["variants_per_window"] = 2
    PIPELINE_CONFIG["max_variants_to_screen"] = 10
    PIPELINE_CONFIG["target_crops"] = ["wheat"]
    print("Configuration updated for testing mode with reduced parameters")

if __name__ == "__main__":
    import json
    config = get_config()
    print("Current Configuration:")
    print(json.dumps(config, indent=2))