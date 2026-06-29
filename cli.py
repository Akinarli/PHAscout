"""
PHAscout Komut Satiri Arayuzu (CLI)
=====================================
Terminal uzerinden hizli ve otomatik analiz yapmak icin kullanilir.

Kullanim:
    python cli.py --accession GCF_000009285.1
    python cli.py --fasta proteins.faa --out report.json
"""

import sys
import click
import logging
from phascout.pipeline import PHAscoutPipeline

# Sadece hatalari veya onemli bilgileri goster
logging.basicConfig(level=logging.WARNING, format="%(message)s")
logger = logging.getLogger("phascout")
logger.setLevel(logging.INFO)


@click.command()
@click.option('--accession', '-a', help='NCBI Assembly Accession (GCF_ / GCA_)')
@click.option('--fasta', '-f', help='Yerel FASTA protein dosyasi')
@click.option('--out', '-o', help='Raporun kaydedilecegi JSON dosyasi')
@click.option('--quiet', '-q', is_flag=True, help='Sadece JSON ciktisi uret (loglari gizle)')
def main(accession, fasta, out, quiet):
    """PHAscout - PHA Uretici Bakteri Tarama Araci"""
    if not accession and not fasta:
        click.echo("HATA: Lutfen --accession veya --fasta belirtin.")
        click.echo("Yardim icin: python cli.py --help")
        sys.exit(1)

    if quiet:
        logger.setLevel(logging.ERROR)

    pipeline = PHAscoutPipeline()

    try:
        report = pipeline.run(accession=accession, fasta_file=fasta)
        
        if out:
            pipeline.reporter.to_json(report, out)
            if not quiet:
                click.echo(f"\n[+] Rapor basariyla kaydedildi: {out}")
        
        if not quiet:
            click.echo("\n" + pipeline.reporter.to_text(report))
        else:
            # Sadece JSON bas (Unix pipelining icin)
            import json
            click.echo(json.dumps(report, ensure_ascii=False))

    except Exception as e:
        click.echo(f"HATA: {e}", err=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
