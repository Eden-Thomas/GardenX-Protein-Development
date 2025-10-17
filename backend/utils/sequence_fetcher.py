"""
sequence_fetcher.py
Fetch real psbA (D1 protein) sequences from UniProt/NCBI
WITH classification labels for contrastive learning
"""

import os
import time
from typing import Dict, Optional, List
from pathlib import Path
import json

# Load environment variables
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

try:
    from Bio import Entrez, SeqIO
    from Bio.Seq import Seq
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("⚠️  Biopython not available. Install with: pip install biopython")


class SequenceFetcher:
    """Fetch protein sequences from public databases with classification labels"""
    
    def __init__(self, email: str = "thomas@edenaglabs.com", ncbi_api_key: Optional[str] = None):
        """
        Args:
            email: Your email (required by NCBI)
            ncbi_api_key: Optional NCBI API key for faster access
        """
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("Biopython is required. Install with: pip install biopython")
        
        # Set up NCBI Entrez
        Entrez.email = email
        if ncbi_api_key:
            Entrez.api_key = ncbi_api_key
            print(f"✅ Using NCBI API key (10 req/sec)")
        else:
            print(f"ℹ️  No NCBI API key (3 req/sec limit)")
        
        self.cache_dir = Path("../data/sequences/cache")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
    
    def fetch_from_uniprot(self, accession: str, organism: str = "") -> Optional[Dict]:
        """Fetch sequence from UniProt by accession ID"""
        import urllib.request
        import urllib.error
        
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        
        try:
            print(f"   Fetching {accession} from UniProt...", end=" ")
            
            with urllib.request.urlopen(url) as response:
                fasta_data = response.read().decode('utf-8')
            
            lines = fasta_data.strip().split('\n')
            header = lines[0]
            sequence = ''.join(lines[1:])
            
            print(f"✅ ({len(sequence)} aa)")
            
            return {
                'accession': accession,
                'organism': organism,
                'sequence': sequence,
                'source': 'UniProt',
                'header': header
            }
            
        except urllib.error.HTTPError as e:
            if e.code == 404:
                print(f"❌ Not found")
            else:
                print(f"❌ Error {e.code}")
            return None
        except Exception as e:
            print(f"❌ Error: {e}")
            return None
    
    def fetch_from_ncbi(self, accession: str, organism: str = "", db: str = "protein") -> Optional[Dict]:
        """Fetch sequence from NCBI by accession ID"""
        
        try:
            print(f"   Fetching {accession} from NCBI...", end=" ")
            
            handle = Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            
            time.sleep(0.34)
            
            print(f"✅ ({len(record.seq)} aa)")
            
            return {
                'accession': accession,
                'organism': organism,
                'sequence': str(record.seq),
                'source': 'NCBI',
                'header': record.description
            }
            
        except Exception as e:
            print(f"❌ Error: {e}")
            return None
    
    def search_ncbi_by_organism(self, organism: str, gene: str = "psbA") -> Optional[str]:
        """Search NCBI for a gene in a specific organism"""
        
        try:
            print(f"   Searching NCBI for {gene} in {organism}...", end=" ")
            
            query = f"{organism}[Organism] AND {gene}[Gene Name] AND chloroplast[filter]"
            
            handle = Entrez.esearch(db="protein", term=query, retmax=5)
            record = Entrez.read(handle)
            handle.close()
            
            time.sleep(0.34)
            
            if record["IdList"]:
                protein_id = record["IdList"][0]
                
                handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
                fasta_record = SeqIO.read(handle, "fasta")
                handle.close()
                
                accession = fasta_record.id.split('|')[0] if '|' in fasta_record.id else fasta_record.id
                
                print(f"✅ Found: {accession}")
                return accession
            else:
                print(f"❌ Not found")
                return None
                
        except Exception as e:
            print(f"❌ Error: {e}")
            return None
    
    def fetch_all_psba_sequences(self) -> Dict[str, Dict]:
        """
        Fetch all psbA sequences with classification labels
        
        Labels:
        - target: Template to improve (maize)
        - positive: Heat-adapted (learn FROM these)
        - negative: Cold-adapted (avoid these traits)
        - neutral: Reference
        """
        
        print("\n" + "="*70)
        print("🔬 FETCHING REAL psbA (D1 PROTEIN) SEQUENCES")
        print("   WITH CONTRASTIVE LEARNING LABELS")
        print("="*70)
        
        sequences = {}
        
        # CORRECTED accessions for complete sequences
        uniprot_targets = {
            # TARGET
            'maize': ('P48183', 'Zea mays', 'target', 28.0),
            
            # POSITIVE (heat-adapted)
            'sorghum': ('P48184', 'Sorghum bicolor', 'positive', 32.0),
            'sugarcane': ('P12212', 'Saccharum officinarum', 'positive', 30.0),  # FIXED
            
            # NEGATIVE (cold-adapted)  
            'wheat': ('P0C433', 'Triticum aestivum', 'negative', 20.0),  # FIXED
            'rice': ('P0C434', 'Oryza sativa', 'negative', 28.0),
            'arabidopsis': ('P56757', 'Arabidopsis thaliana', 'negative', 22.0),  # FIXED
        }
        
        ncbi_targets = {
            # POSITIVE (thermophiles)
            'thermosynechococcus': ('WP_011059484.1', 'Thermosynechococcus elongatus BP-1', 'positive', 55.0),  # FIXED
            'synechococcus': ('WP_011128308.1', 'Synechococcus sp. PCC 7002', 'positive', 38.0),  # FIXED
        }
        
        search_targets = {
            # POSITIVE (heat-adapted C4)
            'setaria': ('Setaria italica', 'positive', 30.0),
            'panicum': ('Panicum miliaceum', 'positive', 30.0),
        }
        
        print("\n📊 Step 1/3: Fetching from UniProt...")
        for species, (accession, organism, label, temp) in uniprot_targets.items():
            result = self.fetch_from_uniprot(accession, organism)
            if result:
                result['label'] = label
                result['temperature'] = temp
                sequences[species] = result
        
        print(f"\n📊 Step 2/3: Fetching from NCBI...")
        for species, (accession, organism, label, temp) in ncbi_targets.items():
            result = self.fetch_from_ncbi(accession, organism)
            if result:
                result['label'] = label
                result['temperature'] = temp
                sequences[species] = result
        
        print(f"\n📊 Step 3/3: Searching NCBI by organism...")
        for species, (organism, label, temp) in search_targets.items():
            accession = self.search_ncbi_by_organism(organism, "psbA")
            if accession:
                result = self.fetch_from_ncbi(accession, organism)
                if result:
                    result['label'] = label
                    result['temperature'] = temp
                    sequences[species] = result
        
        # Validate sequence lengths
        print("\n" + "="*70)
        print(f"✅ FETCHED {len(sequences)} SEQUENCES")
        print("="*70)
        
        valid_sequences = {}
        for species, info in sequences.items():
            length = len(info['sequence'])
            label_emoji = {'target': '🎯', 'positive': '✅', 'negative': '❌', 'neutral': '⚪'}
            emoji = label_emoji.get(info.get('label', 'neutral'), '⚪')
            
            # Validate length (D1 should be ~350-360 aa)
            status = "✅" if 340 <= length <= 370 else "⚠️ "
            
            print(f"{emoji} {species:20s} - {length:3d} aa - {info['source']:8s} - {info['accession']:15s} - {info.get('label', 'neutral'):10s} {status}")
            
            if 340 <= length <= 370:
                valid_sequences[species] = info
            else:
                print(f"   ⚠️  Skipping {species} - unusual length")
        
        print("="*70 + "\n")
        print(f"Valid sequences: {len(valid_sequences)}/{len(sequences)}")
        
        return valid_sequences
    
    def save_sequences(self, sequences: Dict[str, Dict], output_dir: Path = None):
        """Save sequences to files"""
        
        if output_dir is None:
            output_dir = Path("../data/sequences")
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save as FASTA
        fasta_file = output_dir / "psba_sequences_real.fasta"
        with open(fasta_file, 'w') as f:
            for species, info in sequences.items():
                label = info.get('label', 'neutral')
                temp = info.get('temperature', 0)
                f.write(f">{species}|{info['accession']}|{info['organism']}|{label}|{temp}C\n")
                seq = info['sequence']
                for i in range(0, len(seq), 60):
                    f.write(f"{seq[i:i+60]}\n")
        
        print(f"💾 FASTA saved: {fasta_file}")
        
        # Save metadata as JSON
        metadata_file = output_dir / "sequence_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(sequences, f, indent=2)
        
        print(f"💾 Metadata saved: {metadata_file}")
        
        # Save Python dict format for analysis_engine.py
        python_file = output_dir / "sequences_for_analysis_engine.py"
        with open(python_file, 'w') as f:
            f.write("# Real psbA sequences with labels - copy to analysis_engine.py\n\n")
            f.write("REFERENCE_SEQUENCES = {\n")
            for species, info in sequences.items():
                f.write(f"    '{species}': {{\n")
                f.write(f"        'sequence': \"{info['sequence']}\",\n")
                f.write(f"        'temperature': {info.get('temperature', 28.0)},\n")
                f.write(f"        'label': '{info.get('label', 'neutral')}',\n")
                f.write(f"        'photosystem': 'C4',  # Update this\n")
                f.write(f"        'accession': '{info['accession']}',\n")
                f.write(f"        'organism': '{info['organism']}',\n")
                f.write(f"        'description': 'TODO: Add description'\n")
                f.write(f"    }},\n")
            f.write("}\n")
        
        print(f"💾 Python dict saved: {python_file}")
        
        # Generate classification summary
        summary_file = output_dir / "classification_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("SEQUENCE CLASSIFICATION SUMMARY\n")
            f.write("="*70 + "\n\n")
            
            for label in ['target', 'positive', 'negative', 'neutral']:
                label_seqs = [sp for sp, info in sequences.items() if info.get('label') == label]
                if label_seqs:
                    f.write(f"{label.upper()} ({len(label_seqs)} sequences):\n")
                    for sp in label_seqs:
                        f.write(f"  - {sp}: {sequences[sp]['organism']} ({sequences[sp].get('temperature', 0)}°C)\n")
                    f.write("\n")
        
        print(f"💾 Classification summary: {summary_file}")


# Command-line interface
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Fetch psbA sequences from public databases')
    parser.add_argument('--email', type=str, default='thomas@edenaglabs.com',
                       help='Your email (required by NCBI)')
    parser.add_argument('--ncbi-key', type=str, help='NCBI API key (optional)')
    
    args = parser.parse_args()
    
    # Check for API key in environment
    ncbi_key = args.ncbi_key or os.getenv('NCBI_API_KEY')
    
    # Create fetcher
    fetcher = SequenceFetcher(email=args.email, ncbi_api_key=ncbi_key)
    
    # Fetch all sequences
    sequences = fetcher.fetch_all_psba_sequences()
    
    # Save results
    fetcher.save_sequences(sequences)
    
    print("\n✅ SEQUENCE FETCHING COMPLETE!")
    print("\nNext steps:")
    print("1. Review: data/sequences/psba_sequences_real.fasta")
    print("2. Review: data/sequences/classification_summary.txt")
    print("3. Copy sequences from: data/sequences/sequences_for_analysis_engine.py")
    print("4. Update REFERENCE_SEQUENCES in backend/analysis_engine.py")
    print("5. Run Track 1: python3 pipeline/track1_evolutionary.py")