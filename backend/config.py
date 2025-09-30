import os
from pathlib import Path
from dotenv import load_dotenv
import torch

load_dotenv()

class Config:
    # API Keys
    NEUROSNAP_API_KEY = os.getenv("NEUROSNAP_API_KEY", "")
    
    # Compute settings
    COMPUTE_MODE = os.getenv("COMPUTE_MODE", "hybrid")
    
    # Detect available compute device for Apple Silicon
    if torch.cuda.is_available():
        USE_LOCAL_GPU = True
        DEVICE = "cuda"
    elif torch.backends.mps.is_available():
        USE_LOCAL_GPU = True
        DEVICE = "mps"  # Apple Silicon GPU
    else:
        USE_LOCAL_GPU = False
        DEVICE = "cpu"
    
    # Paths
    BASE_DIR = Path(__file__).parent.parent
    DATA_DIR = BASE_DIR / "data"
    SEQUENCES_DIR = DATA_DIR / "sequences"
    RESULTS_DIR = DATA_DIR / "results"
    
    # Model parameters
    THERMOSCORE_WEIGHTS = {
        "esm_fitness": 0.4,
        "composition": 0.3,
        "ogt": 0.3
    }
    
    # Locked positions (1-based)
    LOCKED_POSITIONS = {
        "QB_pocket": [215, 254, 255, 264, 271],
        "PheoD1": [130, 147],
        "TyrZ": [161, 190],
        "OEC_His332": [170, 189, 332, 333, 342, 344]
    }
    
    # Redesign windows for mutations
    REDESIGN_WINDOWS = {
        "window1": range(100, 120),
        "window2": range(200, 220),
        "window3": range(300, 320)
    }
    
    @classmethod
    def create_directories(cls):
        """Create necessary directories if they don't exist"""
        for dir_path in [cls.SEQUENCES_DIR, cls.RESULTS_DIR]:
            if dir_path.exists() and not dir_path.is_dir():
                # If it exists but is not a directory, remove it
                dir_path.unlink()
            dir_path.mkdir(parents=True, exist_ok=True)

# Create directories after class definition
Config.create_directories()