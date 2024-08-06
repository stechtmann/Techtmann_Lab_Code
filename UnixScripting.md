# Help, I want to do something complicated with bash command line!

Basic commands can do a lot on their own (see also other intro to unix documents), but if you want to do something like:
- Repeat the same command across multiple files
- Run a sequence of commands at once
- Modularize and declutter frequently run commands
- Look cool by sending colleagues unintelligible scripts and saying 'run this'
- All of the above at once

*Then you might need to add some loops and/or scripting to your unix toolbelt!*
	
<br />

This tutorial has some brief info and some simple examples for scripting and looping with bash. If you want to learn more (since there is *much* more available to learn), this might also serve as a jumping off point for some googling and stackexchange browsing.

There is also a collection of scripts that might be useful as-is or as examples at the bottom of this page.

## Simple Scripting
Bash scripts (googleable) are essentially the most basic commandline version of an executable (oversimplification). You can put a series of commands in a file that can be run as its own, like so:

```bash
bash yourscript.sh
```
and it will execute everything in that file.

Alternatively, scripts can also be run like this:
```bash
sh yourscript.sh
```
or like this:
```
chmod u+x yourscript.sh

./yourscript.sh
```
---
Let's start by creating a classic hello world script.

<br />

First, we make a new text file using your preferred text editor, such as:
```bash
nano helloworldscript.sh
```
<sub>(traditionally bash scripts use the .sh file extension.)</sub>

<br />

Now that we have the file open, we start with `#!/bin/bash` so that the computer recognizes this is supposed to be a bash script, and follow it with the commands we want executed when the script is run.
```bash
#!/bin/bash

echo "hello world"
```
After the file is saved, our hello world script is ready to go!

<sub>(If you're unfamiliar with nano, you can save+exit by using `ctrl+x` followed by `Y` and `enter` to the subsequent prompts)</sub> 

<br />

Just to be sure, we can test it like so:
```bash
bash helloworldscript.sh
```
which should output
```
hello world
```
<br />

Now, if we want to go crazy with it and actually leverage the fact this is a script, we can add a few extra lines to the script:
```bash
nano helloworldscript.sh
```
```bash
#!/bin/bash

echo "h"
echo "e"
echo "l"
echo "l"
echo "o"
echo ""
echo "w"
echo "o"
echo "r"
echo "l"
echo "d"
```
After saving that, we can test it:
```bash
bash helloworldscript.sh
```
```
h
e
l
l
o

w
o
r
l
d
```
Now we get vertical screen clutter, all with a single command; how useful!

## Variables
Just like other programming languages, you can assign variables in bash (either inside or outside a script).

These will be very helpful for the later parts of this tutorial.

Just to see it in action, let's add some variables to the hello world script:
```bash
nano helloworldscript.sh
```
```bash
#!/bin/bash

#this is a sneaky comment showing comment syntax

#here we define the variable
#note that bash variable names aren't allowed to contain spaces or start with numbers
#variables can be simple text, the output of subcommands, or additional input to the script
TOECHO="variable based hello world"

#here the backticks ` indicate that we want this command run and then stored in the variable
TOECHOTOO=`echo "echoed hello world"`

echo "hello world"
#variables are called with a leading $  (and best practice is to also put them in {})
echo $TOECHO
echo ${TOECHOTOO}
```
and to test:
```
bash helloworldscript.sh
```
```
hello world
variable based hello world
echoed hello world
```

## Looping
<sub>AKA why you're really here</sub>

Similar to other programming languages you can use for and while loops in bash command line, which can be very powerful tools.

It is possible to run loops right on the command line, but it is far easier to manage the syntax if you keep it in a script.

---
#### For loops
As is tradition, a for loop iterates through each element of a collection.

A simple use-case is to make a for loop that goes through each file in the current working directory (which would be wherever you run the script from)
```
nano way_too_verbose_ls.sh
```
```bash
#!/bin/bash

#the for loop syntax creates a variable (here called file) for each element in the list following 'in'
#we then use the variable in the code block beginning with 'do' and ending with 'done'
for file in *
do
		#you can put whatever command(s) you want here
        echo "I have found a file called ${file}"
done
```
and to test:
```
bash way_too_verbose_ls.sh
```
<sub>(I think I'll stick to just using ls)</sub>

---
#### While loops
As is tradition, a while loop repeats as long as a condition remains true. It can be a useful way to have something execute a certain number of times or to create an infinite loop. You can also do such things with a for loop if you wanted.  
<sub>(Infinite loops can be escaped with ctrl+c. Just try not to run them as background processes...)</sub>

```
nano countin.sh
```
```bash
#!/bin/bash
#setting a flag variable
i=1
# -le is bash language for less than or equal to
while [[ $i -le 5 ]] ;
do
        echo $i
        i=$((i + 1))
done
# bash arithmetic syntax is pretty particular
```
testing should result in outputting
```
1
2
3
4
5
```
---
#### More conditionals
If you want to go crazy with it, bash also supports if statements and cases, examples of which may or may not end up here, but can be certainly be found with a bit of googling.

<br />

## Reading inputs
To make your scripts extra modular, you might want to have them read in some extra input (just like a real program!).

There are a couple ways to do this:

#### Command line arguments
Additional arguments will be automatically stored as variables named $1, $2, $3... and so on for each space-separated argument. What you do with these variables is up to you.

```
nano unnecessaryecho.sh
```
```bash
#!/bin/bash
#usage:  `bash unnecessaryecho.sh argToBeEchoed`
echo ${1}
```

#### User input
You can also have the script wait for additional input, and then use that in your script.
For example:
```
nano politeecho.sh
```
```bash
#!/bin/bash

echo "what would you like echoed?"
read TOECHO
echo "Coming right up!"
sleep 3
echo $TOECHO
```

#### Reading a file
Your script can also work from a file.
The recommended way to read in a file (provided here as the first argument), line by line, is like so:
```bash
#!/bin/bash
INPUTFILE=$1

while read line
do
  echo $line
done < "$INPUTFILE"
```
## Useful examples

#### 'Help I want to iterate through each of my sequence files'

```bash
#!/bin/bash

# Assuming sequence files are formatted like samplename_R1_001.fastq.gz and samplename_R1_001.fastq.gz
for file in *_R1*
do
        #This creates a variable from trimming everything before the first 'R1'. Hopefully your files don't have 'R1' in the name.
        PRFX="${file%R1*}"
        #having the prefix separate might be useful downstream
        
        #echoing the filenames just to prove they were found
        #The echo commands can be replaced by something actually useful
        echo "${PRFX}R1_001.fastq.gz"
        echo "${PRFX}R2_001.fastq.gz"
done
```
---
#### 'Help I want to run fastp on all of my sequences'
```bash
#!/bin/bash
# Assuming sequence files are formatted like samplename_R1_001.fastq.gz and samplename_R1_001.fastq.gz
for file in *_R1*
do
        #This creates a variable from trimming everything before the first 'R1'. Hopefully your files don't have 'R1' in the name.
        PRFX="${file%R*}"
        
        #alter fastp flags and arguments to suit your purposes
        fastp --detect_adapter_for_pe --correction --html "fastpreports/${PRFX}.fastp.html" --json "fastpreports/${PRFX}.fastp.json" -i "${PRFX}R1_001.fastq.gz" -I "${PRFX}R2_001.fastq.gz" -m --merged_out "trimmerge/${PRFX}merged.fastq.gz" --include_unmerged --thread 16 --out1 "trimmerge/${PRFX}r1unmerge.fastq.gz" --out2 "trimmerge/${PRFX}r2unmerge.fastq.gz"
done
```
---
#### 'Help I want to grep a list of genes from a series of prokka outputs'
**[WIP]**
```bash
#!/bin/bash
#run this script as 'bash script.sh listofgenes.txt'

#change this to be a path that makes sense (might be complicated)
$PATHTOBINS="."
$GENEFILE="$1"

for genome in $PATHTOBINS*.faa
do
  #this is the recommended read-a-file-line-by-line loop
  while read -r line;
  do
    #modify grep as appropriate
    grep $line $genome
  done < "$GENEFILE"
done
#you might also want to consider some way to aggregate results, since this just spits them out to the command line
```
---
#### 'Help I made a useful bash script and I want it to live here'
Send it to isbigcra@mtu.edu or smtechtm@mtu.edu

