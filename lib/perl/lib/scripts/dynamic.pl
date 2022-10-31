#!/usr/bin/env perl
use Env;
use strict;
use FileHandle;
use Cwd;
use File::Path;
use Sys::Hostname;
use File::Temp qw/ tempfile tempdir /;
no warnings;

my $QUIET="2>/dev/null";
my $VERBOSE=$ENV{VERBOSE_4_DYNAMIC};
my $LOG=$ENV{LOG_4_DYNAMIC};
my $logfile = file2abs("dynamic.log", "new");
our $EXIT_FAILURE=1;
our $EXIT_SUCCESS=0;
our $LAST_COM="";

my %method;
my $method2use;
my $treeF;
my $tree=$ENV{"child_tree_4_TCOFFEE"};
my $thread=$ENV{"child_thread_4_TCOFFEE"};
my $dynamic=$ENV{dynamic_config_4_TCOFFEE};
my $clean;
my $treeFlag;
my $blastFlag;
my $infile;
my $outfile;
my $flush;
my $do_exit=0;
my ($h1, $h2);
my @tmpL;
my $tmpdir = File::Temp->newdir();
my $stderrF="$tmpdir/stderr";
$QUIET="2>$stderrF";
my $CDIR=getcwd();
my $threadFlag4tc;
my $threadFlag4famsa;
my $threadFlag;
my $tcarg;
my $level;

my $QUIET_ENV=$ENV{QUIET_ENV};

if ($QUIET_ENV==1){$QUIET="";}

# -------------------------------------------------------
#     0.    READ COMMAND LINE PARAMS
#           and throw correspondin errors if
#           something is off
# -------------------------------------------------------


for ($a=0; $a<=$#ARGV; $a++)
  {
    if    ($ARGV[$a] eq "-seq"){$infile=file2abs($ARGV[++$a]);}
    elsif ($ARGV[$a] eq "-outfile"){$outfile=file2abs($ARGV[++$a], "new");}
    elsif ($ARGV[$a] eq "-dynamic_config"){
    	$dynamic=file2abs($ARGV[++$a]);
    	if ($VERBOSE){print "\n![dynamic.pl] --- -dynamic_config flag if/else--- $dynamic\n";}
	  }
    elsif ($ARGV[$a] eq "-tree") {$tree=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-method") {$method2use=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-verbose"){$VERBOSE=1; $QUIET="";}
    elsif ($ARGV[$a] eq "-clean"){$clean=1;}
    elsif ($ARGV[$a] eq "-thread"){$thread=$ARGV[++$a]}
    elsif ($ARGV[$a] eq "-tcarg") {$tcarg=file2string($ARGV[++$a]);}
    elsif ($ARGV[$a] eq "-level") {$level=$ARGV[++$a];}
    else
      {
	       add2tcenv($a++,@ARGV);
      }
  }



if ($tree eq "list")
  {
    my $f="$tmpdir/f";
    open (F, ">$f");
    print F ">a\nxxx\n>b\nyyyyy\n";
    close (F);
    print STDOUT ("**** Supported Guide tree modes:\n");
    my_system ("t_coffee -other_pg seq_reformat -in $f -action +seq2dnd list ");
    $do_exit=1;
  }
if ($method2use eq "list")
  {
    my %ml;
    my $listfile="$tmpdir/list";

    $ml{tcoffee}=1;
    $ml{psicoffee}=1;
    $ml{accurate}=1;
    $ml{'3dcoffee'}=1;
    $ml{expresso}=1;
    $ml{clustalo}=1;
    $ml{mafft}=1;
    $ml{famsa}=1;
    $ml{probcons}=1;
    $ml{ginsi}=1;

    print STDOUT ("**** Supported MSA mode:\n");
    my_system ("t_coffee 2>/dev/null | grep _msa > $listfile");
    open (F, $listfile);
    while (<F>)
      {
      	my $l=$_;
      	$l=~/(.*_msa)\s+(.*)/;
      	my $m=$1;
      	my $i="$2\n";
      	if ($m=~/mafftsparsescore/)
    	  {
    	   printf STDOUT "%-20s DOES NOT Support [-tree] -- $i", $m;
    	  }
      	elsif ($m=~/tcoffee/){;}
      	elsif ($m=~/mafft/){;}
      	elsif (!$ml{$m})
    	  {
    	    printf STDOUT "%-20s DOES     Support [-tree] -- $i", $m;
    	  }
      }
    $do_exit=1;
  }
if ($do_exit){my_exit ($CDIR,$EXIT_SUCCESS);}
my $stri=file2string($infile);
my $NSEQ=file2nseq($infile);


if($LOG){
  open (F, $logfile);
  while (<F>){
    print("NSEQ:$NSEQ\n");
    #print("$stri\n");
  }
  close (F);
}


if ($NSEQ==0)
  {
    print "ERROR - No sequences provided [FATAL:dynamic.pl]\n";
    exit ($EXIT_FAILURE);
  }
if (!$outfile)
  {
    ($h1,$outfile)=tempfile();
    push (@tmpL,$outfile);
    $flush=1;
  }


my $master_msa;
if (!($method2use=~/dynamic/)){;}
else
  {
    if (-e $dynamic)
      {
       if ($VERBOSE){print "\n![dynamic.pl] --- -dynamic_config FILE: \n";}
       # Parse dynamic config file
        my @dynamicFile;
        my $index = 0;
      	open (F, $dynamic);
      	while (<F>)
      	  {
      	    my $f=$_;
      	    if ($VERBOSE){print "\n![dynamic.pl] --- FILE content: $f\n";}
            @dynamicFile = split ' ', $f;
            (my $max_nseq = $dynamicFile[1]) =~ s/\s//g;
            if ($VERBOSE){print "\n![dynamic.pl] --- -dynamic_config --- $dynamicFile[0] :: $dynamicFile[1]\n";}

            # Here is the MASTER sequences bucket
            if($index == 0 ){
              print("********* master bucket: $dynamicFile[0] \n");
              $master_msa = $dynamicFile[0];
            }
            # Last case, the one to use for really big buckets
            elsif($max_nseq eq ""){
              print("********* last bucket: $dynamicFile[0] \n");
              $method{$dynamicFile[0]} = "inf";
            }
            # We store the real numbers from the config for the middle cases
            else{
              $method{$dynamicFile[0]} = $dynamicFile[1];
            }
            $index = $index +1;
      	  }
	      close(F);
      }
    else
      {   # default
        	$method{"psicoffee_msa"}=-1;
        	$method{"famsa_msa"}=1000000000;
      }
    # ---------------------------
    #    Select the MSA to use
    # ---------------------------

    # For the maste seqences we use the first
    # method listed
    if ($level==0){
      print("**************************************************************************************************************************************************   MASTER SEQUENCES\n");
      $method2use=$master_msa;
    }else{
      foreach my $name (sort { $method{$a} <=> $method{$b} } keys %method)
        {


            if ($NSEQ<=$method{$name})
              {
                $method2use=$name;
                last;
              }
        }
    }
    print("LEVEL:$level\n");
    print("METHOD:$method2use\n");

    if($LOG){
      open (F, $logfile);
      while (<F>){
        print("LEVEL:$level\n");
        print("METHOD:$method2use\n");
      }
      close (F);
    }

  }


if ($tree)
  {
    ($h2,$treeF)=tempfile();
    my ($h2,$tmptree)=tempfile();
    push (@tmpL,$treeF);
    if ( $tree eq "default"){$treeF=0;}
    elsif ( -e $tree)
      {
	       my_system ("cp $tree $tmptree");
      }
    elsif ($tree eq "master" || $tree eq "main" || $tree eq "parent")
      {
	       if ($ENV{CHILD_TREEF_4_TCOFFEE} && -e $ENV{CHILD_TREEF_4_TCOFFEE})
	        {
      	    my $ctree=$ENV{CHILD_TREEF_4_TCOFFEE};
      	    my_system ("mv $ctree $tmptree");
      	  }
      	else
      	  {
      	    my $master_tree=$ENV{CHILD_TREE_4_TCOFFEE};
      	    my_system ("t_coffee -other_pg seq_reformat -in $master_tree -in2 $infile -action +prune_tree -output newick > $tmptree");
      	  }
      }
    else
      {
	       my_system ("t_coffee -other_pg seq_reformat -in $infile -action +seq2dnd $tree -output newick> $tmptree");
      }

    if ($method2use=~/mafft/)
      {
      	#print "cp $tmptree /Users/cnotredame/.Trash/$$.tmptree\n";
      	#system ("cp $tmptree /Users/cnotredame/.Trash/$$.tmptree");
	      my_system ("t_coffee -other_pg seq_reformat -in $tmptree -output mafftdndmatrix> $treeF");
      }
    else
      {
	      my_system ("mv $tmptree $treeF");
      }
  }
chdir ($tmpdir);

#Collect T-Coffee Command Line
my $CL4tc=get_cl4tc();#will collect from env every CLTCOFEE env variable

# Prepare different tree flag for the different methods
if (!$treeF || $NSEQ<=2){$treeFlag="";}
elsif ( $method2use=~/coffee/ || $method2use=~/accurate/){$treeFlag="-usetree $treeF ";}
elsif ( $method2use=~/clustalo/){$treeFlag="--guidetree-in=$treeF ";}
elsif ( $method2use=~/mafftsparsecore/){;}
elsif ( $method2use=~/mafft/){$treeFlag="--treein $treeF ";}
elsif ( $method2use=~/mafftginsi/){$treeFlag="--treein $treeF ";}
elsif ( $method2use=~/famsa/){$treeFlag="-gt import $treeF ";}
$CL4tc.=" $treeFlag ";

# Prepare different thread flag for different methods
$threadFlag=($thread)?"--thread $thread ":"--thread 1 ";
$threadFlag4tc=($thread)?"-thread $thread ":"-thread 1 ";
$threadFlag4famsa=($thread)?"-t $thread ":"-t 1 ";
$CL4tc.=" $threadFlag4tc ";

# Print Tcoffee ENV
if ($VERBOSE){print "\n![dynamic.pl] --- CL4tc == $CL4tc\n";}


# Launch the command using the prepared flags
my $cmethod=$method2use;
$cmethod=~s/_pair/_msa/;
$cmethod=~s/_msa//;

if ($VERBOSE){print "\n![dynamic.pl] --- cmethod == $cmethod\n";}

if ($cmethod eq "tcoffee"|| $cmethod eq "t_coffee" )
  {
    my_system("cp $CDIR/template_list.txt .");
    my_system("cp $CDIR/*.pdb .");
    my_system ("t_coffee -seq $infile -outfile $outfile -output fasta_aln $CL4tc>/dev/null  $QUIET");
  }
elsif ($cmethod=~/(.*coffee)/ || $cmethod=~/(accurate)/ || $cmethod=~/(expresso)/)
  {
    my $mode=$1;
    my_system ("t_coffee  -mode $mode -seq $infile -outfile $outfile -output fasta_aln $CL4tc >/dev/null  $QUIET");
  }
elsif ($cmethod eq "clustalo")
  {
    my_system ("clustalo -i $infile $treeFlag -o $outfile  --force $threadFlag $QUIET");
    }
elsif ($cmethod eq "mafftginsi")
  {
        my_system ("t_coffee -other_pg seq_reformat -in $treeFlag -input newick -in2 $infile -input2 fasta_seq -action +newick2mafftnewick >> file.mafftnewick");
        print "\n![dynamic.pl][--------------MAFFTGINSI TESTING -------] t_coffee -other_pg seq_reformat -in $treeFlag -input newick -in2 $infile -input2 fasta_seq -action +newick2mafftnewick >> file.mafftnewick \n";

        my_system ("newick2mafft.rb 1.0 file.mafftnewick > file.mafftbinary");
        print "\n \n![dynamic.pl][--------------MAFFTGINSI TESTING -------]newick2mafft.rb 1.0 file.mafftnewick > file.mafftbinary \n";

        my_system ("ginsi --treein file.mafftbinary $infile > $outfile");

  }
elsif ($cmethod =~/sparsecore/)
  {
    my_system ("mafft-sparsecore.rb -i $infile > $outfile $QUIET");
  }
elsif (($cmethod =~/mafft/))
  {
    my $mm;
    my $retree;

    if ( $cmethod eq "mafft" || $cmethod=~/\-/ )
      {
	       $mm=$cmethod;
      }
    elsif (($cmethod=~/mafft(.*)/))
      {
	       $mm="mafft-".$1;
      }
    if ($mm =~/1/)
      {
	       $mm=~s/1/i/;
	       $retree="--retree 1 "
      };
    my_system ("$mm --anysymbol $threadFlag $treeFlag $retree $infile > $outfile $QUIET");
  }
elsif ($method2use=~/famsa/)
  {
    print "\n![dynamic.pl] --- FAMSA DEFAULT\n";
    my_system ("famsa $treeFlag $threadFlag4famsa $infile $outfile >/dev/null $QUIET");
  }
elsif ($method2use=~/famsaUpgma/)
  {
    print "\n![dynamic.pl] --- FAMSA Upgma\n";
    print "\n![dynamic.pl] --- Command: famsa -gt upgma $treeFlag $threadFlag4famsa $infile $outfile >/dev/null $QUIET\n";
    my_system ("famsa -gt upgma $treeFlag $threadFlag4famsa $infile $outfile >/dev/null $QUIET");
  }

elsif ($method2use eq "probcons")
    {
      print "\n![dynamic.pl] --- Command: probcons $infile >  $outfile $QUIET\n";
      my_system ("probcons $infile >  $outfile $QUIET");
    }
else
  {
    if ($treeF)
      {
	       printf (STDERR "WARNING: Method $method2use CANNOT use pre-sepecified guide tree [dynamic.pl]\n");
      }
    my_system ("t_coffee -in $infile -method $method2use -outfile $outfile -output fasta_aln $tcarg -quiet $QUIET");
  }

# Flush output if none provided
if ( ! -e $outfile)
  {
    print "ERROR - No MSA computed - $LAST_COM -- [FATAL:dynamic.pl]\n";
    my_exit ($CDIR,$EXIT_FAILURE);
  }
elsif ( $flush)
 {
   open (F, $outfile);
   while (<F>){print $_;}
   close (F);
 }

# Clean empty files
foreach my $f (@tmpL){unlink($f);}

# Log printing
if ($VERBOSE!=-1)
  {
    open (F, "$stderrF");
    while (<F>)
      {
	      my $l=$_;
	      if ( $VERBOSE || $l=~/WARNING/ || $l=~/ERROR/ || $l=~/INFORMATION/){print stderr "$l";}
      }
    close (F);
  }

my_exit ($CDIR,$EXIT_SUCCESS);



# -------------------------------------------------------
#           HELPER FUNCTIONS
# -------------------------------------------------------

sub file2nseq
  {
    my ($f)=@_;
    my $n=`grep -c '>' $f`;

    return $n;
  }
sub file2abs
    {
      my ($f, $mode)=@_;

      if (!$f || $f=~/^\//){return $f;}
      elsif (!-e $f && $mode eq "new"){return "$CDIR/$f";}
      elsif (!-e $f){return $f;}

      return "$CDIR/$f";
    }
sub file2string
    {
      my ($f)=@_;
      my $s;

      open (F, $f) || return 0;
      while (<F>)
	{
	  $s.=$_;
	}
      close (F);
      chomp($s);
      return $s;
    }

sub get_psicl
      {
	my ($psitrim, $psitrim_mode, $pisN);
	my $cl;

	if ($ENV{psitrim_tree_4_TCOFFEE}){$cl.=" -psitrim_tree=".$ENV{psitrim_tree_4_TCOFFEE}." ";}
	if ($ENV{psitrim_mode_4_TCOFFEE}){$cl.=" -psitrim_mode=".$ENV{psitrim_mode_4_TCOFFEE}." ";}
	if ($ENV{psitrim_4_TCOFFEE}){$cl.=" -psitrim=".$ENV{psitrim_4_TCOFFEE}." ";}
	if ($ENV{psiJ_4_TCOFFEE}){$cl.=" -psiJ=".$ENV{psiJ_4_TCOFFEE}." ";}


	return $cl;
      }

sub get_cl4tc
	{
	  my $cl;

	  foreach my $arg (keys(%ENV))
	    {
	      if ($arg=~/(.*)_4_CLTCOFFEE/)
		{
		  my $name=$1;
		  my $val=$ENV{$arg};
		  if (-e $val){$val=file2abs($val);}


		  if ($val eq "FLAGSET"){$val="";}
		  $cl.="-$name $val ";
		}
	    }
	  if ($VERBOSE){print "\n![dynamic.pl] --- get_psicl --- $cl\n";}
	  return $cl;
	}

sub add2tcenv
	    {
	      my ($p, @argv)=@_;

	      my $flag=$argv[$p];
	      $flag =~s/^-//;
	      my $val =file2abs($argv[$p+1]);
	      my $envv="$flag\_4_CLTCOFFEE";
	      $ENV{$envv}=$val;
	    }

sub my_exit
    {
      my ($dir,$ec)=@_;
      my $a;
      if ($VERBOSE)
	{
	  print "\n![dynamic.pl] --- CDIR: $CDIR\n";
	  print "\n![dynamic.pl] --- Processed $NSEQ\n";
	  print "\n![dynamic.pl] --- ";
	  foreach my $arg (@ARGV)
	    {
	      print "$arg ";
	    }

	  print "\n![dynamic.pl] --- EXIT: $ec ($EXIT_SUCCESS:success, $EXIT_FAILURE:failure)-- Verbose mode -- unset VERBOSE_4_DYNAMIC to turn verbose mode off\n";
	}
      chdir ($dir);
      exit ($ec);
    }

sub my_system
  {
    my ($com)=@_;
    $LAST_COM=$com;

    if ($VERBOSE){print "\n![dynamic.pl] --- SysCall --- $com\n";}

    system ($com);
  }
