"""
download_datasets.py
Automatically download comprehensive D1 datasets from public databases
RUN THIS ONCE to populate your training data
"""

import requests
import time
from Bio import Entrez, SeqIO
from pathlib import Path
from io import StringIO
import json

# CONFIGURE: Set your email (required by NCBI)
Entrez.email = "thomas@gardenx.com"  # Change to your email

class DatasetDownloader:
    """Download D1 sequences from NCBI and UniProt"""
    
    def __init__(self):
        self.data_dir = Path("data/external_databases")
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.sequences = {}
    
    def download_all(self):
        """Run complete download"""
        print("="*70)
        print("DOWNLOADING COMPREHENSIVE D1 DATASETS")
        print("="*70)
        
        # Download thermophile sequences
        print("\n[1/2] Downloading thermophile D1 sequences from NCBI...")
        thermo_seqs = self.download_ncbi_thermophiles()
        
        # Download UniProt sequences
        print("\n[2/2] Downloading curated D1 sequences from UniProt...")
        uniprot_seqs = self.download_uniprot_d1()
        
        # Merge all sequences
        self.sequences.update(thermo_seqs)
        self.sequences.update(uniprot_seqs)
        
        # Save to file
        self.save_sequences()
        
        print("\n" + "="*70)
        print(f"DOWNLOAD COMPLETE! Total sequences: {len(self.sequences)}")
        print("="*70)
        
        return self.sequences
    
    def download_ncbi_thermophiles(self):
        """Download thermophile D1 sequences from NCBI"""
        
        sequences = {}
        
        # Search terms for thermophilic organisms
        organisms = [
            "Thermus thermophilus",
            "Thermococcus kodakarensis",
            "Sulfolobus acidocaldarius",
            "Pyrococcus furiosus",
            "Thermotoga maritima"
        ]
        
        for organism in organisms:
            try:
                print(f"  Searching {organism}...")
                
                # Search
                handle = Entrez.esearch(
                    db="protein",
                    term=f"psbA[Gene] AND {organism}[Organism]",
                    retmax=5
                )
                results = Entrez.read(handle)
                handle.close()
                
                # Download sequences
                for seq_id in results['IdList']:
                    handle = Entrez.efetch(
                        db="protein",
                        id=seq_id,
                        rettype="fasta",
                        retmode="text"
                    )
                    fasta_text = handle.read()
                    handle.close()
                    
                    # Parse
                    record = SeqIO.read(StringIO(fasta_text), "fasta")
                    
                    # Create clean key
                    org_key = organism.lower().replace(' ', '_')
                    sequences[org_key] = str(record.seq)
                    
                    print(f"    ✓ Downloaded {org_key}")
                    
                    time.sleep(0.4)  # Rate limiting
                    
            except Exception as e:
                print(f"    ✗ Error: {e}")
        
        print(f"  Total thermophile sequences: {len(sequences)}")
        return sequences
    
    def download_uniprot_d1(self):
        """Download D1 sequences from UniProt"""
        
        sequences = {}
        
        # UniProt REST API
        query = "gene:psbA AND reviewed:true"
        url = f"https://rest.uniprot.org/uniprotkb/stream?query={query}&format=fasta&size=50"
        
        try:
            print("  Querying UniProt...")
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                # Parse FASTA
                for record in SeqIO.parse(StringIO(response.text), "fasta"):
                    # Extract organism
                    desc = record.description
                    if 'OS=' in desc:
                        organism = desc.split('OS=')[1].split('OX=')[0].strip()
                        
                        # Focus on monocots and key organisms
                        if any(keyword in organism.lower() for keyword in 
                               ['zea', 'triticum', 'oryza', 'sorghum', 'hordeum', 
                                'glycine', 'arabidopsis', 'synechocystis']):
                            
                            org_key = organism.lower().replace(' ', '_')[:30]
                            sequences[org_key] = str(record.seq)
                            print(f"    ✓ {organism}")
                
                print(f"  Total UniProt sequences: {len(sequences)}")
                return sequences
            else:
                print(f"  ✗ UniProt error: {response.status_code}")
                return {}
                
        except Exception as e:
            print(f"  ✗ Error: {e}")
            return {}
    
    def save_sequences(self):
        """Save sequences to Python file for import"""
        
        output_file = self.data_dir / "downloaded_sequences.py"
        
        with open(output_file, 'w') as f:
            f.write('"""\nDownloaded D1 sequences from NCBI and UniProt\n')
            f.write(f'Total sequences: {len(self.sequences)}\n')
            f.write('Auto-generated by download_datasets.py\n"""\n\n')
            
            f.write('ADDITIONAL_SEQUENCES = {\n')
            for org_key, sequence in self.sequences.items():
                # Escape quotes in sequence
                seq_clean = sequence.replace('"', '\\"')
                f.write(f'    "{org_key}": "{seq_clean}",\n')
            f.write('}\n')
        
        # Also save as JSON for backup
        json_file = self.data_dir / "downloaded_sequences.json"
        with open(json_file, 'w') as f:
            json.dump(self.sequences, f, indent=2)
        
        print(f"\n✓ Saved to: {output_file}")
        print(f"✓ Backup: {json_file}")


def main():
    """Run dataset download"""
    
    downloader = DatasetDownloader()
    sequences = downloader.download_all()
    
    print("\n" + "="*70)
    print("NEXT STEPS:")
    print("="*70)
    print("1. Check data/external_databases/downloaded_sequences.py")
    print("2. Run: python backend/integrate_sequences.py")
    print("3. Then run your first analysis!")
    print("="*70)
    
    return sequences


if __name__ == "__main__":
    main()