#%%
# import requests
from abc import abstractclassmethod, abstractmethod
import requests
from collections.abc import Mapping
######################## CLASSE SEQ###############
class Seq :
  def __init__(self, id, description, sequence):
    self._id= id
    self._description=description
    self._sequence=sequence
  @property
  def sequence(self):
    return self._sequence
  @property
  def id (self):
    return self._id
  @property
  def description(self):
    return self._description
    
  def len(self):
    return len(self._sequence)
  def get(self,pos):
   return self._sequence[pos]
  def kmer(self,k):
    seq=self._sequence
    n= len(seq)- k +1
    return set(seq[i:i+k] for i in range(n)) # generateur expression 
    
  def kmers_dist(self, k, seq):
    return (self.kmer(k) & seq.kmer(k))
  def __str__(self):
      return self.id+' ('+self.description+') : '+self.sequence[:10]+'...'
  def __getitem__(self,i):
    return self.get(i) # voir le get



 
#######################CLASSE INVALIDE###############
class InvalidSequenceError(Exception):
  pass
  
#######################CLASSE STANDARCODE ###############
class StandarCode :
  AAs  = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
  Starts = '---M------**--*----M---------------M----------------------------'
  Base1  = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
  Base2  = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
  Base3  = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
  def __init__(self):
    self.condons={}
    n=len(StandarCode.AAs)
    for i in range (n) :
      codon=StandarCode.Base1[i]+StandarCode.Base2[i]+StandarCode.Base3[i]
      aa=StandarCode.AAs[i]
      self.condons[codon]=aa
  def aa(self,codon):
    return self.condons.get(codon.upper(),'*')
    
######################## CLASSE TYPEDSEQ ####################
class TypedSeq(Seq):
  alphabet =None
  def __init__(self,id, description, sequence): # on redef ici pour avoir la vérification comme ça si seq change avec set seq beh ca vérifie quand meme 
    if not self.valid_sequence(self.alphabet, sequence):
      raise InvalidSequenceError('Invalid sequence error')
    super().__init__(id,description,sequence)
  @staticmethod
  def valid_sequence(alphabet,sequence):
    alphabet=set(alphabet.upper()) # tout en maj pour pas de probleme 
    for c in sequence.upper():
      if c not in alphabet :
        return False
    return True
  @classmethod
  def from_ensembl(cls,id):
    record=do_request(ensembl_server, 'sequence/id', id)
    desc=record['desc']
    seq=record['seq']
    #id=record['id']
    molecule=record['molecule']
    return cls(id,desc,seq)
  @staticmethod
  def recode(seq,code):
    csequence = []
    for c in seq:
      csequence.append(code[c])
    return ''.join(csequence)
  def homologs_from_ensembl(self,min_perc_id=0,protein=False):
    ids=ensembl_homology_ids(self=id,min_perc_id=min_perc_id,protein=protein)
    sequences=Sequences()
    sequences.from_ensembl(ids)
    return sequences
    
#######################CLASSE DNASEQ ############
  
class DNASeq (TypedSeq):
  alphabet='actgn' 
  compl={'a':'t','t':'a','c':'g','g':'c','A':'T','T':'A','C':'G','G':'C','n':'n','N':'N'}
  trans={'a':'u','t':'a','c':'g','g':'c','A':'U','T':'A','C':'G','G':'C','n':'n','N':'N'}
 
  def gc (self):
    seq=self._sequence.upper()
    return (seq.count('C')+seq.count('G')) / self.len()
  
  def translate(self):
    code=StandarCode()
    seq=self._sequence.upper()
    protein=''
    n=self.len()-self.len()%3
    for i in range (0, n, 3):
      codon=seq[i:i+3]
      aa=code.aa(codon)
      protein += aa
    return ProteinSeq(self._id, self._description, protein)
  
  
  def reverse_complement(self):
    seq=list(self.complement()._sequence)
    seq=seq.reverse()
    seq=''.join(seq)
    return DNASeq(self._id, self._description, seq)
  def complement(self):
    sequence = self.recode(self._sequence,self.compl)
    return DNASeq(self._id, self._description, sequence)
   
  def transcribe(self):
    sequence=self.recode(self._sequence,self.trans)
    return RNASeq(self._id, self._description, sequence)
######################### CLASSE PROTEINSEQ ##########
class ProteinSeq(Seq):
  alphabet='acdefghiklmnpqrstvwy*'
  
#######################################################
ensembl_server = 'https://rest.ensembl.org'
def do_request(server, service, *args, **kwargs):
    url_params = ''
    for a in args:
        if a is not None:
            url_params += '/' + a
    req = requests.get('%s/%s%s' % (server, service, url_params),
                       params=kwargs,
                       headers={'Content-Type': 'application/json'})
    if not req.ok:
        req.raise_for_status()
    result=req.json()
    return result
def ensembl_homology_ids(id, min_perc_id=0, protein=False):
  hom_response = do_request(ensembl_server, 'homology/id', id)
  homologies = hom_response['data'][0]['homologies']
  id_key = 'id'
  if protein:
    id_key = 'protein_id'
  return [elt['target'][id_key] for elt in homologies
    if elt['source']['perc_id'] > min_perc_id]

class RNASeq (TypedSeq):
  alphabet ='aucgn*'
  btrans={'a':'t','u':'a','c':'g','g':'c','A':'T','U':'A','C':'G','G':'C','n':'n','N':'N'}
  def back_transcribe(self):
    sequence=self.recode(self._sequence,self.btrans)
    return DNASeq(self._id, self._description, sequence)

class Sequences (Mapping) :
  def __init__(self):
    self._sequences={}
    
  def __getitem__(self,key):
    return self._sequences[key]
  def __iter__(self):
    return self._sequences.values().__iter__() # iter(self._sequences.values())

  def __len__(self):
    return self._sequences.__len__() # ou len(self.sequences)
  def add(self,seq):
    if len(self)==0:
      self._first_sequence =seq
    else :
      if type(seq)!=type(self._first_sequence):
        raise InvalidSequenceError("mauvaise seq type")
    self._sequences[seq.id]=seq
  def write_fasta (self,fd):
    for seq in self:
      fd.write('>'+seq.id+'|'+seq.description+'\n')
      fd.write(seq.sequence+'\n\n')
  """def from_fasta (self,fd):
    sequences=Sequences()
    string=fd.read()
    strseq=string.split('>')
    for sseq in strseq:
      if len(sseq)>0 :
        lines=sseq.split('\n')
        id= lines [0].split('|')[0]
        description=lines[0].split('|')[1]
        sequence=''.join(lines[1:])
        seq=Seq(id,description,sequence)
        sequences.add(seq)
    return sequences"""
  @classmethod
  def from_fasta (cls,fd,seq_cls=Seq):
    sequences=cls()
    string=fd.read()
    strseq=string.split('>')
    for sseq in strseq:
      if len(sseq)>0 :
        lines=sseq.split('\n')
        id= lines [0].split('|')[0]
        description=lines[0].split('|')[1]
        sequence=''.join(lines[1:])
        seq=seq_cls(id,description,sequence)
        sequences.add(seq)
    return sequences
  def test_roder(self):
      for seq in self :
        print(seq.id)
  
  def from_ensembl(self,ids):
    
    if len(self)>0:
      seq_cls=type(self._first_sequence)
    else : 
      record=do_request(ensembl_server, 'sequence/id', ids[0])
      if record['molecule']=='protein':
        seq_cls=ProteinSeq
      elif record['molecule']=='dna':
        seq_cls=DNASeq
      elif record['molecule']=='rna':
        seq_cls=RNASeq
    for id in ids:
      seq=seq_cls.from_ensembl(id)
      self.add(seq)
  def kmer_distances(self,k):
    pass

id ='NSG00000157764'
seq=DNASeq.from_ensembl(id)
sequences=seq.homologs_from_ensembl(min_perc_id=99)
sequences.test_order()
import sys
sys.exit(0)
ids= ensembl_homology_ids('ENSG00000157764',min_perc_id=95)
print(ids)
sequences= Sequences()
sequences.from_ensembl(ids)
seq1=DNASeq('id1','desc','ATGCATGCTACGATCG')
seq2=DNASeq('id2','desc','ATGCTAGCTACGTACG')
sequences.add(seq1)
sequences.add(seq2)
seq=sequences['ENSG00000157764']
print(type(seq))

"""for seq in sequences:
  print(seq)"""
with open('seq.fasta','w') as fd :
  sequences.write_fasta(fd)

with open('seq.fasta',"r")as fd: 
  sequences=Sequences.from_fasta(fd,Seq)
  for seq in sequences:
    print(seq)
print(sequences['id2'])
sequences.test_roder()

import sys
sys.exit(0)
  


seq2=DNASeq('idf','description','atgatcgatcgatcgatatcatgcatgcatcgctagcatcgatcgat')
seq3=ProteinSeq('idf','description','atcgtcgagtcgatgcta')
id='ENSG00000157764'
seq=DNASeq.from_ensembl(id)
import sys
seq1=Seq('','','ATCTACGGACTAGCTA')
#seq2=seq
#for k in range(1,seq1.len()):
#  print('k=',k,':',seq1.kmers_dist(k,seq2))
#seq2.kmer(3);
#print(seq.complement()._sequence[:100])
#print(seq.transcribe().back_transcribe()._sequence[:100])
print(seq2)

# %%
