package FASTASequence;

use 5.008;
use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use FASTASequence ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ();
our @EXPORT_OK = ();
our @EXPORT = qw();

our $VERSION = '0.01';


# Preloaded methods go here.
sub new{
  my ($class, $string_text) = @_;
  my ($description,$accession_nr);
  my %f_db_nr = ();
  unless($string_text){
    $string_text = "";
  }
  my $self = {};
  bless($self,$class);

  $string_text =~ s/^(\r?\n)+//;
  $string_text =~ s/^>+/>/;
  my ($description_line,$sequence) = split(/\n/,$string_text,2);
  $sequence =~ s/\r?\n//g;
  unless($description_line =~ /^>/){
    $self->{accession_nr} = "";
  }
  else{
    if($description_line =~ /^>gi\|/){
      my ($gi,$number,$db);
      $description_line =~ s/^>//;
      ($gi,$number,$db,$accession_nr,$description) = split(/\|/,$description_line);
    }
    elsif($description_line =~ /^>sp\|/){
      $description_line =~ s/^>//;
      $description = (split(/\s/,$description_line,2))[1];
      my $desc = (split(/\s/,$description_line,2))[0];
      $accession_nr = (split(/\|/,$description_line))[1];
    }
    elsif($description_line =~ /^>[XY\d+]/){
      $description_line =~ s/>//;
      chomp $description_line;
      $description = (split(/\s/,$description_line,3))[-1];
      $accession_nr = (split(/\s/,$description_line,3))[0];
    }
    elsif($description_line =~ /^>[0-9A-Z_]+\s?/){
      $description_line =~ s/^>//;
      #-------------------------------------------------#
      # IPI-Sequences                                   #
      #-------------------------------------------------#
      if($description_line =~ /^IPI:/){
        # split only at first whitespace and take first element
        my $foreign_numbers = (split(/\s/,$description_line,2))[0];
        $description = (split(/\s/,$description_line,3))[2];
        my @foreign_acs = split(/\|/,$foreign_numbers);
        # cross-references to other databases
        foreach my $f_ac(@foreign_acs){
          my ($key, $value) = split(/:/,$f_ac);
          $f_db_nr{$key} = $value;
        }
	unless($f_db_nr{'SWISS-PROT'}){
	  $f_db_nr{'SWISS-PROT'} = "NULL";
	}
	unless($f_db_nr{'ENSEMBL'}){
	  $f_db_nr{'ENSEMBL'} = "NULL";
	}
	unless($f_db_nr{'REFSEQ_XP'}){
	  $f_db_nr{'REFSEQ_XP'} = "NULL";
	}
	unless($f_db_nr{'TREMBL'}){
	  $f_db_nr{'TREMBL'} = "NULL";
	}
        $accession_nr = $f_db_nr{'IPI'};
	delete $f_db_nr{IPI};
      }
      #-----------------------------------------#
      # format begins with accession-nr         #
      #-----------------------------------------#
      elsif($description_line =~ /^[A-Z][0-9][A-Z0-9]{3}[0-9][\s\|]/){
        $description_line =~ s/^>//;
        if($description_line =~ /\|/){
          ($accession_nr, $description) = split(/\|/,$description_line,2);
        }
	else{
          ($accession_nr, $description) = split(/\s/,$description_line,2);
	}
      }
      elsif($description_line =~ /^[A-Z0-9_]+$/){
        $description_line =~ s/^>//;
        chomp $description_line;
	$accession_nr = $description_line;
      }
      else{
        $description_line =~ s/^>//;
        my $seq_id;
        ($seq_id, $accession_nr,$description)= split(/\s/,$description_line,3);
      }
    }
  }

  $accession_nr =~ s/^>//;
  $accession_nr =~ s/\.\d//;
  $sequence =~ s/[^A-Z]//g;

  $self->{text}         = $sequence;
  $self->{accession_nr} = $accession_nr;
  $self->{description}  = $description;
  $self->{seq_length}   = length($sequence);
  $self->{dbrefs}       = \%f_db_nr;
  $self->{crc64}        = $self->_crc64();

  return $self;
}#end new

sub getSequence{
  my ($class) = @_;
  return $class->{text};
}# end getText

sub getSequenceLength{
  my ($class) = @_;
  return $class->{seq_length};
}# end getSequenceLength

sub getAccessionNr{
  my ($class) = @_;
  return $class->{accession_nr};
}# end of getAccessionNr

sub getDescription{
  my ($class) = @_;
  return $class->{description}
}# end getDescription

sub getCrc64{
  my ($class) = @_;
  return $class->{crc64};
}# end getCrc64

sub getDBRefs{
  my ($class) = @_;
  return $class->{dbrefs};
}# end getDBRefs

sub allIndexesOf{
  my ($self,$search) = @_;
  my $i = 1;
  my $index = 0;
  my @indices = ();
  while($i != -1){
    $index = index($self->{text},$search,$index);
    push(@indices,$index) unless ($index == -1);
    $i = $index;
    $index++;
  }
  return \@indices;
}# end allIndicesOf

sub _crc64 {
  my ($self)     = @_;
  my $text = $self->{text};
  use constant EXP => 0xd8000000;
  my @highCrcTable = 256;
  my @lowCrcTable  = 256;
  my $initialized  = ();
  my $low          = 0;
  my $high         = 0;

  unless($initialized) {
    $initialized = 1;
    for my $i(0..255) {
      my $low_part  = $i;
      my $high_part = 0;
      for my $j(0..7) {
        my $flag = $low_part & 1; # rflag ist für alle ungeraden zahlen 1
        $low_part >>= 1;# um ein bit nach rechts verschieben
        $low_part |= (1 << 31) if $high_part & 1; # bitweises oder mit 2147483648 (), wenn $parth ungerade
        $high_part >>= 1; # um ein bit nach rechtsverschieben
        $high_part ^= EXP if $flag;
      }
      $highCrcTable[$i] = $high_part;
      $lowCrcTable[$i]  = $low_part;
    }
  }

  foreach (split '', $text) {
    my $shr = ($high & 0xFF) << 24;
    my $tmph = $high >> 8;
    my $tmpl = ($low >> 8) | $shr;
    my $index = ($low ^ (unpack "C", $_)) & 0xFF;
    $high = $tmph ^ $highCrcTable[$index];
    $low  = $tmpl ^ $lowCrcTable[$index];
  }
  return sprintf("%08X%08X", $high, $low);
}# end crc64

sub seq2file{
  my ($self,$file,$args_ref) = @_;
  open(W_SEQUENCE,">$file") or die "Can't open $file: $!\n";
  print W_SEQUENCE ">",$self->{accession_nr};
  foreach my $dbkey(keys(%{$self->{dbrefs}})){
  print W_SEQUENCE "|".$dbkey.":".$self->{dbrefs}->{$dbkey};
  }
  print W_SEQUENCE " ",$self->{description},"\n";
  print W_SEQUENCE $self->{text},"\n";
  close W_SEQUENCE;
}# end seq2file

sub getFASTA{
  my ($self)     = @_;
  my $fasta = ">".$self->{accession_nr};
  foreach my $dbkey(keys(%{$self->{dbrefs}})){
  $fasta .= "|".$dbkey.":".$self->{dbrefs}->{$dbkey};
  }
  $fasta .= " ".$self->{description}."\n";
  $fasta .= $self->{text}."\n";
  return $fasta;
}# end getFASTA

sub addDBRef{
  my ($self,$db,$dbref) = @_;
  if($self->{dbrefs}->{$db} && ($self->{dbrefs}->{$db} ne 'NULL')){
    $self->{dbrefs}->{$db} .= ";".$dbref;
  }
  else{
    $self->{dbrefs}->{$db} = $dbref;
  }
}# end addDBRef

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

FASTASequence - Perl extension Biooinformatics

=head1 SYNOPSIS

  use FASTASequence;
  my $fasta = qq~>sp|P01815|HV2B_HUMAN Ig heavy chain V-II region COR - Homo sapiens (Human).
QVTLRESGPALVKPTQTLTLTCTFSGFSLSSTGMCVGWIRQPPGKGLEWLARIDWDDDKY
YNTSLETRLTISKDTSRNQVVLTMDPVDTATYYCARITVIPAPAGYMDVWGRGTPVTVSS
  ~;
  my $seq = FASTASequence->new($fasta);

=head1 ABSTRACT

  This should be the abstract for FASTASequence.
  The abstract is used when making PPD (Perl Package Description) files.
  If you don't want an ABSTRACT you should also edit Makefile.PL to
  remove the ABSTRACT_FROM option.

=head1 DESCRIPTION

This perl module is a simple utility to simplify the job of bioinformatics.
It parses several information about a given FASTA-Sequence such as:

=over 10

=item * accession number

=item * description

=item * sequence itself

=item * length of sequence

=item * crc64 checksum (as it is used by SWISS-PROT)

=back

=head2 METHODS

=head3 new

=head3 getAccessionNr

	my $accession = $seq->getAccessionNr();

returns the AccessionNr of the FASTA-Sequence

=head3 getDescription

	my $description = $seq->getDescription();

returns the description standing in the first line of the
FASTA-format (without the accession number)

=head3 getSequence

	my $sequence = $seq->getSequence();

returns the sequence

=head3 getCrc64

	my $crc64_checksum = $seq->getCrc64();

returns the crc64 checksum of the sequence. This checksum
corresponds with the crc64 checksum of SWISS-PROT

=head3 addDBRef

	$seq->addDBRef(DB, REFERENCE_AC);

DB is the name of the referenced database

REFERENCE_AC is the accession number in the referenced database

=head3 seq2file

	$seq->seq2file(FILENAME, OPTIONS);

FILENAME is the path of the file where the sequence has to be stored.

OPTIONS is a hash, which contains the options:

=over 10

=item -o

overwrite the output-file if the file already exists. true / false
Default: true

=back

=head3 allIndexesOf

	my $indexes = $seq->allIndexesOf(EXPR);

returns a reference on an array, which contains all indexes of
EXPR in the sequence

=head3 getSequenceLength

	my $length = $seq->getSequenceLength();

returns the length of the sequence

=head3 getDBRefs

	my $hashref = $seq->getDBRefs();

returns a hashreference. The hash contains all references
	hashref = {'SWISS-PROT' => 'P01815'},

=head3 getFASTA

	my $fasta_sequence = $seq->getFASTA();

returns the sequence in FASTA-format

=head2 EXAMPLE

	use FASTASequence;
	my $fasta = qq~>sp|P01815|HV2B_HUMAN Ig heavy chain V-II region COR - Homo sapiens (Human).
	QVTLRESGPALVKPTQTLTLTCTFSGFSLSSTGMCVGWIRQPPGKGLEWLARIDWDDDKY
	YNTSLETRLTISKDTSRNQVVLTMDPVDTATYYCARITVIPAPAGYMDVWGRGTPVTVSS
	~;

	my $seq = FASTASequence->new($fasta);

	print 'The sequence of '.$seq->getAccessionNr().' is '.$seq->getSequence(),"\n";
	print 'This sequence contains '.scalar($seq->allIndexesOf('C').' times Cystein at the following positions:';
	print $_+1.', ' for(@{$seq->allIndexesOf('C')});

=head1 ADDITIONAL INFORMATION

=head3 accepted formats

This module can parse the following formats:

=over 4

=item >P02656 APC3_HUMAN Apolipoprotein C-III precursor (Apo-CIII).

=item >IPI:IPI00166553|REFSEQ_XP:XP_290586|ENSEMBL:ENSP00000331094|TREMBL:Q8N3H0 T Hypothetical protein

=item >sp|P01815|HV2B_HUMAN Ig heavy chain V-II region COR - Homo sapiens (Human).

=back

=head3 structure

The structure of the hash for the example is:

	$VAR1 = {
	         'seq_length' => 120,
	         'accession_nr' => 'P01815',
	         'text' => 'QVTLRESGPALVKPTQTLTLTCTFSGFSLSSTGMCVGWIRQPPGKGLEWLARIDWDDDKYYNTSLETRLTISKDTSRNQVVLTMDPVDTATYYCARITVIPAPAGYMDVWGRGTPVTVSS',
	         'crc64' => '158A8B29AE7EEB98',
	         'dbrefs' => {},
	         'description' => 'Ig heavy chain V-II region COR - Homo sapiens (Human).'
	       }

if you miss something please contact me.

=head1 BUGS

There is no bug known. If you experienced any problems, please contact me.

=head1 SEE ALSO

http://perl-modules.renee-baecker.de

the crc64-routine is based on the
SWISS::CRC64
module.

=head1 AUTHOR

Renee Baecker, E<lt>module@renee-baecker.deE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2004 by Renee Baecker

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. 

=cut
