from qiime2.plugin import Metadata

HEALTHY_SPECIES_DEFAULT = [
    "s__Alistipes_senegalensis",
    "s__Bacteroidales_bacterium_ph8",
    "s__Bifidobacterium_adolescentis",
    "s__Bifidobacterium_angulatum",
    "s__Bifidobacterium_catenulatum",
    "s__Lachnospiraceae_bacterium_8_1_57FAA",
    "s__Sutterella_wadsworthensis"]

NON_HEALTHY_SPECIES_DEFAULT = [
    "s__Anaerotruncus_colihominis",
    "s__Atopobium_parvulum",
    "s__Bifidobacterium_dentium",
    "s__Blautia_producta",
    "s__candidate_division_TM7_single_cell_isolate_TM7c",
    "s__Clostridiales_bacterium_1_7_47FAA",
    "s__Clostridium_asparagiforme",
    "s__Clostridium_bolteae",
    "s__Clostridium_citroniae",
    "s__Clostridium_clostridioforme",
    "s__Clostridium_hathewayi",
    "s__Clostridium_nexile",
    "s__Clostridium_ramosum",
    "s__Clostridium_symbiosum",
    "s__Eggerthella_lenta",
    "s__Erysipelotrichaceae_bacterium_2_2_44A",
    "s__Flavonifractor_plautii",
    "s__Fusobacterium_nucleatum",
    "s__Gemella_morbillorum",
    "s__Gemella_sanguinis",
    "s__Granulicatella_adiacens",
    "s__Holdemania_filiformis",
    "s__Klebsiella_pneumoniae",
    "s__Lachnospiraceae_bacterium_1_4_56FAA",
    "s__Lachnospiraceae_bacterium_2_1_58FAA",
    "s__Lachnospiraceae_bacterium_3_1_57FAA_CT1",
    "s__Lachnospiraceae_bacterium_5_1_57FAA",
    "s__Lachnospiraceae_bacterium_9_1_43BFAA",
    "s__Lactobacillus_salivarius",
    "s__Peptostreptococcus_stomatis",
    "s__Ruminococcaceae_bacterium_D16",
    "s__Ruminococcus_gnavus",
    "s__Solobacterium_moorei",
    "s__Streptococcus_anginosus",
    "s__Streptococcus_australis",
    "s__Streptococcus_gordonii",
    "s__Streptococcus_infantis",
    "s__Streptococcus_mitis_oralis_pneumoniae",
    "s__Streptococcus_sanguinis",
    "s__Streptococcus_vestibularis",
    "s__Subdoligranulum_sp_4_3_54A2FAA",
    "s__Subdoligranulum_variabile",
    "s__Veillonella_atypica"]


def _load_and_validate_species(list_healthy: str = None,
                               list_non_healthy: str = None):

    healthy_species = _load_file(list_healthy) \
        if list_healthy else HEALTHY_SPECIES_DEFAULT
    non_healthy_species = _load_file(list_non_healthy) \
        if list_non_healthy else NON_HEALTHY_SPECIES_DEFAULT

    if not healthy_species:
        raise ValueError('Healthy species list is empty!')
    if not non_healthy_species:
        raise ValueError('Non-healthy species list is empty!')

    return healthy_species, non_healthy_species


def _load_file(file: str = None):
    with open(file, 'r') as f:
        return list(map(lambda x: x.strip(), f.readlines()))


def _load_metadata(metadata: Metadata = None):
    if not metadata:
        raise ValueError('Metadata parameter not provided!')
    metadata = metadata.to_dataframe()
    return metadata
