"""
visualizations.py
Generate data for interactive visualizations
"""

import numpy as np
from typing import Dict, List
import json

class ConservationPlotter:
    """Generate conservation plot data"""
    
    @staticmethod
    def generate_plot_data(conservation_scores: np.ndarray, window_size: int = 10) -> Dict:
        """Generate smoothed conservation plot data"""
        # Smooth scores with moving average
        smoothed = np.convolve(conservation_scores, np.ones(window_size)/window_size, mode='valid')
        
        positions = list(range(1, len(conservation_scores) + 1))
        
        return {
            'labels': positions[::5],  # Every 5th position for readability
            'raw_data': conservation_scores.tolist(),
            'smoothed_data': smoothed.tolist(),
            'threshold_high': [0.9] * len(positions),
            'threshold_medium': [0.7] * len(positions)
        }

class CorrelationHeatmap:
    """Generate correlation heatmap data"""
    
    @staticmethod
    def generate_heatmap(sequences: Dict[str, str]) -> Dict:
        """Generate pairwise identity heatmap"""
        species_list = list(sequences.keys())
        n = len(species_list)
        
        # Calculate pairwise identities
        identity_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(n):
                if i == j:
                    identity_matrix[i][j] = 100
                else:
                    seq1 = sequences[species_list[i]]
                    seq2 = sequences[species_list[j]]
                    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
                    identity_matrix[i][j] = (matches / len(seq1)) * 100
        
        return {
            'species': species_list,
            'matrix': identity_matrix.tolist(),
            'colorscale': 'YlOrRd'
        }

class MutationImpactChart:
    """Generate mutation impact visualization data"""
    
    @staticmethod
    def generate_chart_data(recommendations: List[Dict]) -> Dict:
        """Generate bubble chart data for mutation impacts"""
        data = {
            'mutations': [],
            'x_values': [],  # Confidence
            'y_values': [],  # Predicted ΔΔG
            'sizes': [],     # Thermophile prevalence
            'colors': []     # Conservation score
        }
        
        for rec in recommendations[:15]:  # Top 15
            mutation_label = f"{rec['from']}{rec['position']}{rec['to']}"
            data['mutations'].append(mutation_label)
            data['x_values'].append(rec['confidence'])
            data['y_values'].append(abs(rec['delta_stability']))
            data['sizes'].append(rec['thermophile_prevalence'])
            data['colors'].append(rec['conservation_score'])
        
        return data

class FeatureImportancePlot:
    """Generate feature importance plots"""
    
    @staticmethod
    def generate_importance_data(feature_importance: Dict[str, float]) -> Dict:
        """Generate bar chart data for feature importance"""
        # Sort by importance
        sorted_features = sorted(feature_importance.items(), key=lambda x: x[1], reverse=True)
        
        # Take top 10
        top_features = sorted_features[:10]
        
        return {
            'labels': [feat[0].replace('_', ' ').title() for feat in top_features],
            'values': [feat[1] for feat in top_features],
            'colors': ['#667eea'] * len(top_features)
        }

class SequenceLogoGenerator:
    """Generate sequence logo data"""
    
    @staticmethod
    def generate_logo_data(aligned_sequences: Dict[str, str], start: int = 0, end: int = 50) -> Dict:
        """Generate sequence logo data for a region"""
        sequences = list(aligned_sequences.values())
        logo_data = []
        
        for pos in range(start, min(end, len(sequences[0]))):
            column = [seq[pos] for seq in sequences if seq[pos] != '-']
            
            if not column:
                continue
            
            # Calculate frequencies
            aa_counts = {}
            for aa in column:
                aa_counts[aa] = aa_counts.get(aa, 0) + 1
            
            total = len(column)
            frequencies = {aa: count/total for aa, count in aa_counts.items()}
            
            # Calculate information content (bits)
            entropy = -sum(f * np.log2(f) for f in frequencies.values() if f > 0)
            max_entropy = np.log2(20)  # 20 amino acids
            information = max_entropy - entropy
            
            logo_data.append({
                'position': pos + 1,
                'residues': frequencies,
                'information': information
            })
        
        return logo_data