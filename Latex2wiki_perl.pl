# This Perl 5 script translates LaTeX into wiki text. It is definitely not
# complete, but can handle much more than the original Python script,
# in particular macros can be resolved. It can be used as a CGI
# script or on the command line (when certain parts are commented or
# uncommented, see below "BEGINNING of"). The output will be HTML in 
# both cases.

# Written by E Mueller
#
# The GNU license of the original Python script (Latex2wiki.py by
# Maxime Biais, cf. http://wiki.loria.fr/wiki/Latex2wiki, see also
# http://www.kataplop.net/pub/info/projets/latex2twiki) applies to
# this script, too:
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
use strict;
## This part is used within a CGI script - comment it if the script is used on the command line
## BEGINNING of CGI script part
#  use CGI qw/:standard/;
#  autoEscape(undef);
#  my $tex = '';
#  my $documentation = '';
#  if (param())
#  {
#    $tex = param('field_name');
#    $documentation = param('info');
#  }
#  if ($documentation)
#  {
#    print header,
#          start_html('Documentation'),
#          h1('Documentation');
#  }
#  else
#  {
#    print header,
#          start_html('TeX-Wiki Conversion'),
#          h1('TeX to Wiki Conversion'),
#          p,
#          a({-href=>url(-relative=>1)."?info=1"}, "Documentation"),
#          start_form,
#          "TeX or LaTeX text ",p,textarea(-name=>'field_name',
#                           -default=>$tex,
#                           -override=>1,
#                           -columns=>50,
#                           -rows=>20),p,
#          submit,        
#          end_form,
#          hr;
#  }
#  if (param()) 
#  {
#    my $error;
#    my $wiki;
#    $wiki = process($tex, \$error, $documentation);
#    if ($error)
#    {
#      print $error,p;
#    }
#    else
#    {
#      if ($documentation)
#      {
#        print '<pre>', $wiki, '</pre>';
#      }
#      else
#      {
#        print "Paste this Wiki text",p,textarea(-override=>1,-columns=>50, -rows=>20, -default=>$wiki), 
#              hr;
#      }
#    }
#  }
#  print end_html;
## END of CGI script part

## This part is used for command line - comment it if the script is used as a CGI script
## BEGINNING of command line part
my $error;
unless(@ARGV)
{
  print "usage: 'perl $0 -h' for help, 'perl $0 foo.tex' reads file 'foo.tex'";
  exit (0);
}
if ($ARGV[0] eq '-h')
{
  $_ = process("", \$error, 1); # Documentation
}
else
{
  $_ = process(join('', (<>)), \$error, 0);
}
if ($error)
{
  print STDERR $error;
}
else
{
  print "<html><body><pre>";
  print $_;
  print "</pre></body></html>";
}
## END of command line part


# obtain text until the end of the current block or environment:
# parameter 0: respect environments, pos $_ is before end of environment/closing brace
# parameter 1: only respect braces; pos $_ is after closing brace
sub getBrace($)
{
  my $level = 1;
  my $position = pos($_);
  if ($_[0])
  {
    while  (m/\G.*?(?<!\\)([{}])/gcs)
    {
      $level += ($1 eq '{')? 1:-1;
      return substr($_, $position, pos($_)-1-$position) if ($level == 0);
    }
  }
  else
  { 
    while (m/\G.*?(?<!\\)(\{|\}|\\begin|\\end)/gcs)
    {
      $level += ($1 eq '{' or $1 eq '\begin')? 1:-1;
      if ($level == 0)
      {
        pos($_) -= length($1);
        return substr($_, $position, pos($_)-$position);
      }
    }
  }
  return substr($_, $position); # if braces do not match, simply return end of string
}  

# obtain next TeX token
sub getToken()
{
  return $1 if (m/\G([\w\s])/gc); # normal character: one token
  return $1 if (m/\G(\\\w+)\s*/gc); # control sequence
  return $1 if (m/\G(\\[\#\$\&\\\{\}\%])/gc); # escaped control characters
  return getBrace(1) if (m/\G\{/gc); # text in braces: text until closing brace
  return ""; # end
}


# parses through the tokens of the text, replacing macros
# second parameter: verbatim contents
sub replace_macros(\%\@)
{
  my ($macros, $verbatim) = @_;
  pos $_ = 0;
  my $position = 0;
  my $token;
  my $star; # star version of macro, e. g. \hspace*
  my $distance = length $_; # distance from current position to end: for control of infinite recursion
  my $recursion_count = 0;
  my $start_verbatim = 0; # no verbatim before \begin{document}
  while (m/\G.*?((\\\w+)(\*?)\s*)/sgc)
  {
    $token = $2;
    $star = $3;
    $position = pos($_)-length($1); # starting position of control sequence
    if ($token eq '\begin')
    {
      my $type = getToken();
      if ($type eq "document")
      {
        $start_verbatim = 1;
      }
      elsif ($start_verbatim and $type eq "verbatim")
      {
        pos($_) = $position;
        if (s/\G.*?\{verbatim\}(.*?)\\end\s*\{verbatim\}/\x00\x01/sc)
        {
          my $verb_text = $1;
          chomp($verb_text);
          $verb_text =~ s#\n#</nowiki>\n <nowiki>#g;
          push @$verbatim, " <nowiki>$verb_text</nowiki>\n";
          pos($_) = $position+2;
          next;
        }
      }
    }
    elsif ($start_verbatim and $token eq '\verb')
    {
      pos($_) = $position;
      if (s/\G\s*\\verb\\*?\s*(.)(.*?)\1/\x00\x02/sc)
      {
        push @$verbatim, $2;
        pos($_) = $position+2;
        next;
      }
    }
    if (length($_)-$position < $distance) # distance to end of string
    {
      $recursion_count = 0;
      $distance = length($_)-$position;
    }
    else
    {
      $recursion_count ++;
      if ($recursion_count > 1000)
      {
         return "Maybe infinite recursion at macro $token in document";
      }
    }
    my $replacement = "";
    # if only non-star version of macro exists and star appears: star does not belong to token
    if ($star and exists $$macros{$token} and not exists $$macros{$token.$star})
    {
      $star = '';
      pos($_) = $position+length($token); # position directly after token without star - star will be next token
    }
    $token .= $star;
    if (exists $$macros{$token})
    {
      my $paramcount = $$macros{$token}[0];
      $replacement = $$macros{$token}[1];
      if ($paramcount > 0)
      {
        my @parameters = ("");
        # optionaler Parameter?
        if ($#{$$macros{$token}} > 1 and m/\G\[(.*?)\]/gc)
        {
          $parameters[0] = $1;
          $replacement = $$macros{$token}[2];
        }
        while ($paramcount-- > 0)
        {
          push @parameters, getToken();
        }
        $replacement =~ s/\#(\d)/$parameters[$1]/g;
      }
      elsif ($paramcount < 0)
      {
        my $parameter = getBrace(0);
        $replacement =~ s/\#1/$parameter/g;
      }
    }
    elsif ($token eq '\def')
    {
      my $paramcount = 0;
      $token = getToken();
      $paramcount = $1 if (m/\G[\#\d]+(\d)/gc);
      @{$$macros{$token}} = ( $paramcount, getToken() );
    }     
    elsif ($token eq '\newcommand' or $token eq '\renewcommand')
    {
      my $paramcount = 0;
      $token = getToken();
      $paramcount = $1 if (m/\G\[(\d)\]/gc);
      @{$$macros{$token}} = ( $paramcount, getToken() );
    }
    else
    {
      next; # no replacement
    }
    substr($_, $position, (pos $_)-$position) = $replacement;
    pos $_ = $position;
  }
  return ""; # no infinite recursion
}



# itemize/enumeration environment
sub items($$$\$)
{
  if ($_[0] eq "item")
  {
    return "\n".${$_[3]};
  }
  if ($_[1] eq "begin") # in Wiki 
  {
    ${$_[3]} .= ($_[2] eq "itemize") ? '*' : '#';
  }
  else
  {
    chop ${$_[3]};
    return "\n\n" unless (${$_[3]}); # at the end of first enumeration/itemization reset counters using empty line
  }
  return "";
}

# convert content of tabular environment to Wiki table with small borders. It will not consider any further options.
sub tabelle($)
{
  my $inhalt = $_[0];
  $inhalt =~ s#\\\\[^\&]*$##;
  $inhalt =~ s#[\r\n]# #g;
  $inhalt =~ s#[\r\n]# #g;
  $inhalt =~ s#\&#\n\|#g;
  $inhalt =~ s#\\\\#\n\|-\n\|#g;
  return "\n\x00\x03 border='1' cellspacing='0'\n\|$inhalt\n\x00\x04";
}

# actual conversion function.
# parameters: 
# - TeX text
# - reference to error text
# - documentation mode?
# result: Wiki text in HTML format
sub process($\$$)
{
  # standard macros. For each macro number of parameters and replacement text (with parameters
  # as in TeX). If number of parameters is -1, then take everything until the end of the current
  # block or environment. An optional third parameter indicates a replacement in which the
  # that the contents of the optional parameter [..] can be inserted as #0
  my %macros = ( '\centerline' => [1,'<center>#1</center>'],
  '\textit' => [1,"''#1''"],
  '\texttt' => [1,"=#1="],
  '\textbf' => [1,"'\|\|\|\|apos;'#1'\|\|\|\|apos;'"], # the browser condenses two literal apostrophes to a quotation mark
  '\bf' => [-1, "'\|\|\|\|apos;'#1'\|\|\|\|apos;'"],
  '\it' => [-1, "'\|\|\|\|apos;#1\|\|\|\|apos;'"],
  '\emph' => [1,"'\|\|\|\|apos;#1\|\|\|\|apos;'"],
  '\title' => [1,"= #1 ="],
  '\ref' => [1, "---#1---"],
  '\ss' => [0, "||||szlig;"],
  '\aa' => [0, "||||aring;"],
  '\AA' => [0, "||||Aring;"],
  '\ae' => [0, "||||aelig;"],
  '\AE' => [0, "||||AElig;"],
  '\o' => [0, "||||oslash;"],
  '\O' => [0, "||||Oslash;"],
  '\dq' => [0, "||||quot;"],
  '\c' => [1, "||||#1cedil;"],
  '\H' => [1, "||||#1uml;"],
  '\settowidth' => [2, ""],
  '\addtolength' => [2, ""],
  '\setlength' => [2, ""],
  '\vspace' => [1, ""],
  '\hspace' => [1, ""],
  '\vspace*' => [1, ""],
  '\hspace*' => [1, ""],
  '\cite' => [1, "[#1]", "[#1 #0]"]
   ); 
  # The result will be written into this array and finally joined
  my @parts;
  if ($_[2]) # Documentation
  {
    push @parts, (<<'EOF');
This script converts TeX and LaTeX input to wiki text which can be pasted from a web browser. (The literal output is
always HTML in order to avoid encoding problems).
It resolves macros which are defined in the input (e. g. \def\foo{..}, \def\foo#1{..}, \newcommand\foo{...}, 
\renewcommand\foo{...}, \newcommand\foo[1]{...}, \renewcommand\foo[1]{...}), and detects infinite recursions.
(unlike TeX/LaTeX, the definition will be valid until the end of the document, not
just to the end of the enclosing block).
It respects quotes and umlauts from german.sty (if text contains \input german or \usepackage{german}).
Of course, it cannot consider any other referenced files (\input, \usepackage etc.).
It supports umlauts and accents \", \', \`, \^, \~ for e. g. ||||Auml;, ||||Aacute;, ||||Agrave;, ||||Acirc;, ||||Atilde;.
It replaces verbatim, itemize, center, and enumeration environments and \verb and
\section, \subsection, etc. with corresponding Wiki elements. It handles simple tabular environments
(replaces them with a table, but does not support more sophisticated features such as alignment and \hline)
It handles %, \%, \$, \par, ~, \~, \\, \_, --, ---, replaces \@. with . and removes \/ and replaces math mode ($...$, 
$$...$$, \(...\), \[...\]; and equation, align, eqnarray environments) with proper Wiki <math>...</math> context.
It replaces the following macros as if they were defined as (note that the actual outcome may look differently in 
the browser because the browser will replace the &...; commands):
EOF
    foreach my $macro (sort keys %macros)
    {
      next if $macros{$macro}[0] < 0;
      my $replacement = $macros{$macro}[1];
      my $parameter = $macros{$macro}[0] > 0 ? "\[$macros{$macro}[0]\]" : "";
      $replacement =~ s/\|\|\|\|/\&/g;
      push @parts, "\\newcommand ${macro}${parameter}\{$replacement\}", ($macros{$macro}[2] ? " (the optional \[...\] argument is supported)" : "") ,  "\n";
    }
    push @parts, "It handles the following TeX commands, where #1 denotes everything until the end of the current block:\n";
    foreach my $macro (sort keys %macros)
    {
      next if $macros{$macro}[0] >= 0;
      my $replacement = $macros{$macro}[1];
      $replacement =~ s/\|\|\|\|/\&/g;
      push @parts, "${macro} -> $replacement\n";
    }
    push @parts, "\nAll undefined macros and environments in text mode will just be removed, but their parameters will remain.\n",
                 "Note that \\cite will be properly handled only if the parameter is an URL\n",
                 "\nWritten by E M\|\|\|\|uuml;ller in 2009-2010, inspired by some Python scripts";
  }
  else # normal processing
  {
    my @verbatim; # contents of verbatim environment (will be extracted from text and re-inserted in the end)
    $_ = $_[0];  
    $_ .= "\n";
    # delete comments
    s/(?<!\\)\%.*?\n//g;
    # replace percent signs
    s/\\%/\%/g;
    # replace macros
    my $error = replace_macros(%macros, @verbatim);
    if ($error)
    {
      ${$_[1]} = $error;
      return;
    }
    undef %macros;
    # replace multiple empty lines or \par by <br/>
    s#\n\s*\n+\s*|\\par\s*#<br/>\n#g;
    # replace multiple space with one space
    s#\s{2,}# #sg;
    # centered text
    s#\\begin\{center\}(\s*?\n)?#<center>#g;
    s#\\end\{center\}(\s*?\n\s*)?#</center>#g;
    s#\\begin\{tabular\}\{.*?\}(.*)\\end\{tabular\}#&tabelle($1)#sge;
    # german?
    my $german;
    $german = 1 if (/\\input\s+german/ or /\\usepackage\s*\{.*?german.*?\}/);
    
    # items
    my $beginning; # beginning of wiki line with item
    s#\n?\\((item\b)(\[.*?\])?|(begin|end)\s*\{(itemize|enumerate)\})\s*#&items($2, $4, $5, \$beginning)#gse;
  
    # replace \i by i in accents (in TeX \i is used to remove dot above i):
    s#(\\["'`´^])\{?\\([Ii])\}?#$1$2#g;
    # accents and umlauts
    s#\\"\{?([aeiouAEIOU])\}?#\|\|\|\|$1uml;#g;
    s#\\'\{?([aeiouAEIOUyY])\}?#\|\|\|\|$1acute;#g;
    s#\\`\{?([aeiouAEIOU])\}?#\|\|\|\|$1grave;#g;
    s#\\\^\{?([aeiouAEIOU])\}?#\|\|\|\|$1circ;#g;
    s#\\['´]\{?([aeiouAEIOU])\}?#\|\|\|\|$1acute;#g;
    s#\\\~\{?([ANOano])\}?#\|\|\|\|$1tilde;#g;
    # quotes and ampersand
    s#``|''|´´#\|\|\|\|quot;#g;
    s#\\&#\|\|\|\|amp;#g;
    
    # umlauts with german.sty:
    if ($german)
    {
      s#\"([aeiouAEIOU])#\|\|\|\|$1uml;#g;
      s#\"s#\|\|\|\|szlig;#g;
      s#\"['`"]#\|\|\|\|quot;#g;
    } 
    # break text into mathematical and non-mathematical parts. Substitutions only in non-mathematical parts
    # preparation: convert into TeX: $ and $$
    s#\\[\(\)]#\$#g;
    s#\\[\[\]]#\$\$#g;
    s#\\(begin|end)\s*\{equation\*?\}(\s*?\n)?#\$\$#g;
    s#\\begin\{(align|eqnarray)\*?\}(\s*?\n)?#\$\$\\begin{matrix}#g;
    s#\\end\{(align|eqnarray)\*?\}(\s*?\n)?#\\end{matrix}\$\$#g;
    # replace display math by centering and normal math
    my $toggledisplay = 0;
    s#\$\$#$toggledisplay++; ($toggledisplay%2) ? '<center>$': '$</center>'#egs;
    # splitting into mathematical and non-mathematical mode
    @parts = split /(?<!\\)\$/;
    $_ = "";
    # removing beginning of document:
    if ($parts[0] =~ /\\document(class|style)/)
    {
      $parts[0] =~ s/^.*\\begin\{document\}(\s*\n)?//s;
    }
    for (my $i=0; $i < @parts; $i+=2)
    { 
      # section, subsection, ...
      $parts[$i] =~ s#\\(sub)*section\*?\{(.*?)\}#my $level="=" x (2+length($1)/3);"\n$level$2$level\n"#eg;
      # line breaks
      $parts[$i] =~ s#\\\\\n?\s*#<br/>#g;
      # delete remaining environments
      $parts[$i] =~ s#\\begin\s*\{minipage\}\s*(\[\w*\])\s*\{.*?\}\s*##g;
      $parts[$i] =~ s#\\(begin|end)\s*\{.*?\}(\s*?\n)?##sg;
      # replace multiple dashes by one dash
      $parts[$i] =~ s#---?#-#g;
      # replace obligatory space by normal space 
      $parts[$i] =~ s#\\ # #g; 
      # underlines
      $parts[$i] =~ s#\\_#_#g; 
      # dollar signs
      $parts[$i] =~ s#\\\$#\$#g; 
      # remove macros
      $parts[$i] =~ s#\\\w+\*?\s*##sg; 
      # remove braces
      $parts[$i] =~ s#[\{\}]##g;
      # replace unescaped tildes with non-breaking spaces
      $parts[$i] =~ s#^\~#\&nbsp;#g;
      $parts[$i] =~ s#([^\\])\~#$1\&nbsp;#g;
      $parts[$i] =~ s#\\\~#\~#g;
      # replace \@. by .
      $parts[$i] =~ s#\\\@\.#.#g;
      # remove \/
      $parts[$i] =~ s#\\/##g;
      # insert verbatim content
      $parts[$i] =~ s#\x00\x01#"<br/><code>".shift(@verbatim)."</code><br/>"#ge;
      $parts[$i] =~ s#\x00\x02#"<code>".shift(@verbatim)."</code>"#ge;
      $parts[$i] =~ s#\x00\x03#\{\|#g;
      $parts[$i] =~ s#\x00\x04#\|\}#g;
    }
    # remove empty lines in the beginning
    $parts[0] =~ s#^(\s*<br/>)*\s*##sg;
    # remove empty lines in the end
    $parts[-1] =~ s#(\s*<br/>)*\s*$##sg;
    # insert math delimiters
    for (my $i=1; $i < @parts; $i+=2)
    {
      $parts[$i] = "<math>$parts[$i]</math>";
    }
  }
  my $result =join('', @parts);
  undef @parts;
  # replace HTML entities (and replace |||| by &)
  $result =~ s/\&/&amp;/g;
  $result =~ s/\|\|\|\|/\&/g;
  $result =~ s/</&lt;/g;
  $result =~ s/>/&gt;/g;
  return $result;
}

