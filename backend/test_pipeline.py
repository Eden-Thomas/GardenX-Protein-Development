import asyncio
from pipeline.master_pipeline_v3 import MasterPipelineV3Sequential
from config import Config

async def main():
    config = Config()
    pipeline = MasterPipelineV3Sequential(config, client=None)
    
    # Corn D1 sequence
    wt_sequence = 'MTAILERRESESLWGRFCNWITGTENRLYIGWFGVLMIPTLLTATSVFIIAFIAAPPVDIDGIREPVSGSLLYGNNIITGAVVPSSNAILVHFYPIWEAASLDEWLYNGGPYQLIIFHFLLGASCYMGREWELSFRLGMRPWICVAYSAPLASAFAVFLIYPIGQGSFSDGMPLGISGTFNFMIVFQAEHNILMHPFHMLGVAGVFGGSLLFSAMHGSLVTSSLIRETTENESANEGYRFGQEEETYNIVAHAGYGRLLFQYASFNNSRSLHFFLAAWPVVGIWFTALGISTMAFNLNGFNFNQSVVDSQGRVINTWADIINRANLGMEVMHERNAHNFPLDLAAEG'
    
    print('🌽 Engineering thermostable corn D1 protein')
    print(f'   Target: +27°C improvement (28°C → 55°C)')
    print(f'   Method: True sequential Track 1 → Track 2 → Track 3')
    
    results = await pipeline.run(
        wt_sequence=wt_sequence,
        target_temp=55,
        n_variants=5,
        validate_all=False
    )

if __name__ == "__main__":
    asyncio.run(main())