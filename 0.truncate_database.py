import os, sys, django
sys.path.append('/data/django/rgenome_net')
os.environ['DJANGO_SETTINGS_MODULE'] = 'rgenome_net.settings'
from django.conf import settings
django.setup()

from cas_database.models import *

raw_input("This script deletes all records from database. Proceed? (Ctrl+C to cancel)")

print("Deleting Target_Transcripts...")
Target_Transcript.objects.all().delete()
print("Deleting Targets...")
Target.objects.all().delete()
print("Deleting Offtargets...")
Offtarget.objects.all().delete()
print("Deleting Guides...")
Guide.objects.all().delete()
print("Deleting Sequences...")
Sequence.objects.all().delete()
print("Deleting Transcripts...")
Transcript.objects.all().delete()
print("Deleting CDSranges...")
CDSrange.objects.all().delete()
print("Deleting Genes...")
Gene.objects.all().delete()
print("Deleting Chromosomes...")
Chromosome.objects.all().delete()
print("Deleting Organisms...")
Organism.objects.all().delete()
