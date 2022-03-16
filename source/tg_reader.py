import gzip
import sys

#
# accepts fq / fq.gz / fa / fa.gz
#
# handles fa files with newlines in read sequence
#

class TG_Reader:
	def __init__(self, input_filename, verbose=True):
		self.verbose = verbose
		fnl = input_filename.lower()
		#
		if fnl[-3:] == '.fq' or fnl[-6:] == '.fq.gz' or fnl[-6:] == '.fastq' or fnl[-9:] == '.fastq.gz':
			self.filetype = 'FASTQ'
		elif fnl[-3:] == '.fa' or fnl[-6:] == '.fa.gz' or fnl[-6:] == '.fasta' or fnl[-9:] == '.fasta.gz':
			self.filetype = 'FASTA'
		else:
			print('Error: unknown file suffix given to TG_Reader:')
			print(input_filename)
			exit(1)
		#
		if fnl[-3:] == '.gz':
			if self.verbose:
				print('getting reads from gzipped ' + self.filetype + '...')
			self.f = gzip.open(input_filename, 'rt')
		else:
			if self.verbose:
				print('getting reads from ' + self.filetype + '...')
			self.f = open(input_filename, 'r')
		#
		self.buffer = []
		self.current_readname = None

	#
	# returns (readname, readsequence, qualitysequence)
	#
	def get_next_read(self):
		if self.filetype == 'FASTQ':
			my_name = self.f.readline().strip()[1:]
			if not my_name:
				return ('','','')
			my_read = self.f.readline().strip()
			skip    = self.f.readline().strip()
			my_qual = self.f.readline().strip()
			return (my_name, my_read, my_qual)
		#
		elif self.filetype == 'FASTA':
			if self.current_readname == None:
				self.current_readname = self.f.readline().strip()[1:]
			if not self.current_readname:
				return ('','','')
			hit_eof = False
			while True:
				my_dat = self.f.readline().strip()
				if not my_dat:
					hit_eof = True
					break
				self.buffer.append(my_dat)
				if '>' in self.buffer[-1]:
					break
			if hit_eof:
				out_dat = (self.current_readname, ''.join(self.buffer), '')
				self.current_readname = None
				self.buffer = []
			else:
				out_dat = (self.current_readname, ''.join(self.buffer[:-1]), '')
				self.current_readname = self.buffer[-1][1:]
				self.buffer = []
			return out_dat

	def close(self):
		self.f.close()

if __name__ == '__main__':
	#
	IN_READS_TEST = sys.argv[1]
	my_reader = TG_Reader(IN_READS_TEST)
	while True:
		read_dat = my_reader.get_next_read()
		if not read_dat[0]:
			break
		print(read_dat)
	my_reader.close()
