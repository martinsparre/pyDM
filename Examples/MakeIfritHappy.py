import commands

GadgetFiles = commands.getoutput('find . -name \*_\*').split('\n')

for f in GadgetFiles:
    commands.getoutput( 'ln -s ' + f + ' ' + f.replace('_','_0') )
