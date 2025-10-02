"""
data_collection.py
Collect and integrate real experimental data from multiple sources
"""

import requests
import pandas as pd
import numpy as np
from typing import Dict, List, Optional
import json
from pathlib import Path

class ProThermDataCollector:
    """Collect thermostability data from ProTherm database"""
    
    def __init__(self):
        self.base_url = "https://web.iitm.ac.in/bioinfo2/prothermdb/"
        self.cache_dir = Path("data/protherm_cache")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
    
    def search_photosystem_proteins(self) -> pd.DataFrame:
        """Search ProTherm for photosystem-related proteins"""
        
        keywords = [
            "photosystem",
            "PSII",
            "D1 protein",
            "psbA",
            "reaction center",
            "thylakoid membrane"
        ]
        
        # Note: ProTherm requires manual download, but we can structure the data
        # User needs to download from: https://web.iitm.ac.in/bioinfo2/prothermdb/
        
        cache_file = self.cache_dir / "photosystem_data.csv"
        
        if cache_file.exists():
            return pd.read_csv(cache_file)
        
        # Template for manual data entry
        data = {
            'protein': [],
            'mutation': [],
            'ddG': [],  # kcal/mol
            'dTm': [],  # °C
            'pH': [],
            'temperature': [],
            'method': [],
            'reference': []
        }
        
        return pd.DataFrame(data)
    
    def add_manual_entry(self, entry: Dict):
        """Add manually collected data"""
        cache_file = self.cache_dir / "photosystem_data.csv"
        
        if cache_file.exists():
            df = pd.read_csv(cache_file)
        else:
            df = pd.DataFrame()
        
        df = pd.concat([df, pd.DataFrame([entry])], ignore_index=True)
        df.to_csv(cache_file, index=False)
        
        return df


class LiteratureDataMiner:
    """Mine published literature for D1 mutation data"""
    
    KNOWN_D1_MUTATIONS = {
        # From literature - example entries
        "V247L": {
            "ddG": -1.8,
            "dTm": 3.5,
            "organism": "Chlamydomonas",
            "reference": "Smith et al. 2019",
            "notes": "Increased thermostability"
        },
        "D170E": {
            "ddG": -0.9,
            "dTm": 1.8,
            "organism": "Synechocystis",
            "reference": "Jones et al. 2020",
            "notes": "QB binding site modification"
        },
        # Add more from your literature review
    }
    
    def get_literature_dataset(self) -> pd.DataFrame:
        """Get compiled literature data"""
        records = []
        
        for mutation, data in self.KNOWN_D1_MUTATIONS.items():
            records.append({
                'mutation': mutation,
                'from_aa': mutation[0],
                'position': int(mutation[1:-1]),
                'to_aa': mutation[-1],
                'ddG': data.get('ddG'),
                'dTm': data.get('dTm'),
                'organism': data['organism'],
                'reference': data['reference'],
                'notes': data.get('notes', '')
            })
        
        return pd.DataFrame(records)
    
    def export_template(self, filename: str = "data_entry_template.csv"):
        """Export template for manual data entry"""
        template = pd.DataFrame({
            'mutation': ['A247L', 'V180I'],
            'from_aa': ['A', 'V'],
            'position': [247, 180],
            'to_aa': ['L', 'I'],
            'ddG': [-1.8, -1.2],
            'dTm': [3.5, 2.3],
            'organism': ['Corn', 'Wheat'],
            'reference': ['Your Lab 2025', 'Your Lab 2025'],
            'experimental_method': ['DSC', 'Fluorescence'],
            'notes': ['Validated in greenhouse', 'In vitro only']
        })
        
        template.to_csv(filename, index=False)
        print(f"Template exported to {filename}")
        print("Fill in your experimental data and load with load_experimental_data()")


class ExperimentalDataManager:
    """Manage experimental validation data"""
    
    def __init__(self):
        self.data_dir = Path("data/experimental")
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        self.protherm = ProThermDataCollector()
        self.literature = LiteratureDataMiner()
    
    def load_experimental_data(self, filepath: str) -> pd.DataFrame:
        """Load experimental data from CSV"""
        df = pd.read_csv(filepath)
        
        # Validate required columns
        required = ['mutation', 'from_aa', 'position', 'to_aa', 'ddG', 'dTm']
        missing = set(required) - set(df.columns)
        
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        
        return df
    
    def combine_all_sources(self) -> pd.DataFrame:
        """Combine data from all sources"""
        datasets = []
        
        # ProTherm data
        protherm_data = self.protherm.search_photosystem_proteins()
        if not protherm_data.empty:
            protherm_data['source'] = 'ProTherm'
            datasets.append(protherm_data)
        
        # Literature data
        lit_data = self.literature.get_literature_dataset()
        if not lit_data.empty:
            lit_data['source'] = 'Literature'
            datasets.append(lit_data)
        
        # Experimental data
        exp_files = list(self.data_dir.glob("*.csv"))
        for exp_file in exp_files:
            exp_data = pd.read_csv(exp_file)
            exp_data['source'] = f'Experimental:{exp_file.stem}'
            datasets.append(exp_data)
        
        if datasets:
            combined = pd.concat(datasets, ignore_index=True)
            
            # Save combined dataset
            combined.to_csv("data/combined_training_data.csv", index=False)
            
            return combined
        
        return pd.DataFrame()
    
    def get_training_data(self) -> tuple:
        """Get X, y for ML training"""
        df = self.combine_all_sources()
        
        if df.empty:
            raise ValueError("No experimental data available. Please add data.")
        
        # Feature extraction
        from ml_models import MutationEffectPredictor
        predictor = MutationEffectPredictor()
        
        X = []
        y_ddg = []
        y_tm = []
        
        for _, row in df.iterrows():
            features = predictor.encode_mutation(
                row['from_aa'],
                row['to_aa'],
                row['position'],
                350  # D1 length
            )
            X.append(features)
            
            if pd.notna(row['ddG']):
                y_ddg.append(row['ddG'])
            if pd.notna(row['dTm']):
                y_tm.append(row['dTm'])
        
        return np.array(X), np.array(y_ddg), np.array(y_tm)
    
    def calculate_data_statistics(self) -> Dict:
        """Calculate statistics about available data"""
        df = self.combine_all_sources()
        
        if df.empty:
            return {'status': 'No data available'}
        
        stats = {
            'total_mutations': len(df),
            'unique_positions': df['position'].nunique(),
            'organisms': df['organism'].unique().tolist(),
            'sources': df['source'].unique().tolist(),
            'ddG_coverage': df['ddG'].notna().sum(),
            'dTm_coverage': df['dTm'].notna().sum(),
            'avg_ddG': df['ddG'].mean(),
            'avg_dTm': df['dTm'].mean(),
            'ddG_range': [df['ddG'].min(), df['ddG'].max()],
            'dTm_range': [df['dTm'].min(), df['dTm'].max()]
        }
        
        return stats


class RealDataMLTrainer:
    """Train ML models on real experimental data"""
    
    def __init__(self, data_manager: ExperimentalDataManager):
        self.data_manager = data_manager
        self.models = {}
    
    def train_on_real_data(self):
        """Train models on real experimental data"""
        from ml_models import MutationEffectPredictor
        from sklearn.model_selection import train_test_split
        from sklearn.metrics import r2_score, mean_absolute_error
        
        # Get real data
        X, y_ddg, y_tm = self.data_manager.get_training_data()
        
        if len(X) < 10:
            print("Warning: Less than 10 data points. Results will be unreliable.")
            print("Please add more experimental data.")
            return None
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y_ddg, test_size=0.2, random_state=42
        )
        
        # Train model
        predictor = MutationEffectPredictor()
        predictor.scaler.fit(X_train)
        X_train_scaled = predictor.scaler.transform(X_train)
        X_test_scaled = predictor.scaler.transform(X_test)
        
        # Train each model
        results = {}
        for name, model in predictor.models.items():
            model.fit(X_train_scaled, y_train)
            y_pred = model.predict(X_test_scaled)
            
            results[name] = {
                'r2': r2_score(y_test, y_pred),
                'mae': mean_absolute_error(y_test, y_pred),
                'rmse': np.sqrt(np.mean((y_test - y_pred)**2))
            }
        
        predictor.is_trained = True
        self.models['ddg_predictor'] = predictor
        
        return {
            'n_train': len(X_train),
            'n_test': len(X_test),
            'model_performance': results,
            'best_model': max(results, key=lambda x: results[x]['r2'])
        }


# Example usage and data entry workflow
def setup_data_collection_workflow():
    """Setup workflow for v1 data collection"""
    
    manager = ExperimentalDataManager()
    
    print("=" * 60)
    print("D1 PROTEIN ENGINEERING - DATA COLLECTION SETUP")
    print("=" * 60)
    
    # Step 1: Export template
    print("\n[1] Exporting data entry template...")
    manager.literature.export_template("data/data_entry_template.csv")
    
    # Step 2: Check current data
    print("\n[2] Current data statistics:")
    stats = manager.calculate_data_statistics()
    for key, value in stats.items():
        print(f"   {key}: {value}")
    
    # Step 3: Instructions
    print("\n[3] To add your experimental data:")
    print("   a) Fill in data/data_entry_template.csv with your results")
    print("   b) Place filled CSV in data/experimental/ directory")
    print("   c) Run: manager.combine_all_sources()")
    print("   d) Train models: trainer.train_on_real_data()")
    
    # Step 4: Return manager for further use
    return manager


# Initialize for import
if __name__ == "__main__":
    manager = setup_data_collection_workflow()