from Bio import SeqIO


records = SeqIO.to_dict(SeqIO.parse('data/LANL-HIV.fasta', 'fasta'))
extraction_information = [
  { 'input_id': 'B.KR.2005.05CSR3.DQ837381', 'output_id': 'related_1' },
  { 'input_id': 'B.US.2012.505_0122a.WG08.MG196700', 'output_id': 'related_2' },
  { 'input_id': 'B.CH.2002.HIV_CH_BID-V3527_2002.JQ403021', 'output_id': 'diverged_1' },
  { 'input_id': 'P.FR.2006.RBF168.GQ328744', 'output_id': 'diverged_2' }
]

for info in extraction_information:
  record = records[info['input_id']]
  SeqIO.write(record, 'data/%s.fasta' % info['output_id'], 'fasta')

