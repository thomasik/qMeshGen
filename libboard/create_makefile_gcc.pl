#!/usr/bin/perl
##

open(MAKEFILE,">Makefile");

print MAKEFILE "CXX       = g++\n";
print MAKEFILE "DIR       = linux\n";
print MAKEFILE "CXXFLAGS  = -DUNIX_MODE -std=c++11 -pedantic -I. -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -Wno-unknown-pragmas -Wno-char-subscripts\n";
print MAKEFILE "LDFLAGS   = -O2\n\n";
print MAKEFILE "MESHLIB   = libboard.a\n";
print MAKEFILE "SRCS      =";

@srcs = @ARGV;

$i = 0;
foreach $file (@srcs){
    @file_headers = scan_headers($file);
    foreach $header (@file_headers){
        if(!exists $h_headers{$header}){
            @sub_headers = scan_headers($header);
            $h_headers{$header} = join("\t", @sub_headers);
        }
    }
    $c_headers{$file} = join("\t", @file_headers);
    print MAKEFILE " $file";
    if(++$i == 2){
        print MAKEFILE "\\\n\t";
        $i = 0;
    }
}
print MAKEFILE "\n";

$i = 0;
print MAKEFILE "\nOBJS      = ";
foreach $file (@srcs){
    $obj_file = "\$(DIR)/".$file;
    $obj_file =~ s/\.cpp/\.o/;
    $object_files{$file} = $obj_file;
    print MAKEFILE " $obj_file";
    if(++$i == 2){
        print MAKEFILE "\\\n\t";
        $i = 0;
    }
}
print MAKEFILE "\n";

print MAKEFILE 'default: $(DIR)/$(MESHLIB)', "\n\n";
print MAKEFILE '$(DIR)/$(MESHLIB): $(OBJS)';
print MAKEFILE "\n\tar rcs \$(DIR)/\$(MESHLIB) \$(OBJS) \n\n";
print MAKEFILE "clean:\n\trm -f \$(OBJS) \$(DIR)/\$(MESHLIB)\n\n";

# TODO follow headers

foreach $file (@srcs){
    print MAKEFILE "$object_files{$file}: $file";
    @headers = split(/\t/, $c_headers{$file});
    foreach $h (@headers){
        print MAKEFILE " $h";
    }
    print MAKEFILE "\n\t\${CXX} \$(CXXFLAGS) -c $file -o $object_files{$file}\n\n";
}

sub scan_headers($) {
    my($filename) = @_;
    my(@local_headers) = undef;

    print "Scanning $filename...\n";

    open(INFILE,"$filename");
    while($line = <INFILE>){
        if($line =~ /include \"(.*)\"/){
            push(@local_headers, $1) if -e $1;
        }
    }
    print "--- @local_headers ---\n";
    return @local_headers;
}
