# team_with_fetch.py
import requests

team = [
    {
        "name": "Opeoluwa",
        "slack_username": "Eadencre8ives",
        "country": "Nigeria",
        "hobby": "Troubleshooting & Movies",
        "affiliations": ["Microbiology", "Bioinformatics", "Molecular Biology"],
        "favorite_gene": "BRCA1",
        # prefer to fetch; fallback used if fetch fails
        "favorite_gene_accession": "NM_007294.4",  # BRCA1 mRNA RefSeq (example)
        "favorite_gene_sequence": None
    }
]

NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

def fetch_sequence_from_ncbi(accession, rettype="fasta"):
    """Fetch sequence from NCBI using efetch. Returns FASTA string or None."""
    params = {
        "db": "nucleotide",
        "id": accession,
        "rettype": rettype,
        "retmode": "text"
    }
    try:
        resp = requests.get(NCBI_EFETCH_URL, params=params, timeout=20)
        resp.raise_for_status()
        return resp.text
    except Exception as e:
        print(f"[Warning] Could not fetch {accession} from NCBI: {e}")
        return None

def fasta_to_sequence(fasta_text):
    """Convert FASTA text to raw sequence (single-line, uppercase)."""
    if not fasta_text:
        return None
    lines = fasta_text.splitlines()
    seq_lines = [ln.strip() for ln in lines if ln and not ln.startswith(">")]
    return "".join(seq_lines).upper()

def prepare_team_sequences(team_list, fetch=True):
    for member in team_list:
        accession = member.get("favorite_gene_accession")
        if fetch and accession:
            fasta = fetch_sequence_from_ncbi(accession)
            seq = fasta_to_sequence(fasta)
            if seq:
                member["favorite_gene_sequence"] = seq
            else:
                member["favorite_gene_sequence"] = "SEQUENCE_NOT_AVAILABLE"
        else:
            if not member.get("favorite_gene_sequence"):
                member["favorite_gene_sequence"] = "SEQUENCE_NOT_PROVIDED"

def print_team_info(team_list):
    for m in team_list:
        print(f"Name:           {m['name']}")
        print(f"Slack username: {m['slack_username']}")
        print(f"Country:        {m['country']}")
        print(f"Hobby:          {m['hobby']}")
        print(f"Affiliations:   {', '.join(m['affiliations'])}")
        print(f"Favorite gene:  {m['favorite_gene']}")
        print("Gene sequence (first 200 bases):")
        seq = m.get("favorite_gene_sequence") or "N/A"
        # print a snippet for readability; show full if short
        print(seq[:200] + ("..." if len(seq) > 200 else ""))
        print("-" * 60)

if __name__ == "__main__":
    # Attempt to fetch sequences; set fetch=False to skip network call
    prepare_team_sequences(team, fetch=True)
    print_team_info(team)
