import os, sys
from glob import glob
import re
from mutagen.mp3 import MP3
 
#### DEBUG THIS ####
def compare_songs(br1, br2, size1, size2):
   if (br1 < br2):
      return 2
   elif (br1 == br2):
      if (size1 < size2):
	     return 2
   return 1

dir_name = raw_input("Which directory: ")
band_name = raw_input("Enter artist name: ")
print dir_name
print band_name

#CD to band directory
os.chdir(dir_name)

#delete un-needed image files
ext_list = ['gif', 'jpg', 'jpeg','png','txt', 'PNG', 'JPG', 'GIF']
for ext in ext_list: 
   for filename in glob("*.%s" %ext):
      os.remove(filename)

if (os.path.exists("../%s.log"% band_name)):
   os.remove("../%s.log" %band_name)	  

files = glob("*.*")

print "There are %d files to be sorted in this directory" % len(files)
bitrate_dict = {}
name_dict = {}
count_dict = {}

fp = open("../%s.log" %band_name, 'w')

for song in files:
   old_name = song	#it will soon be renamed
   print "checking: ", old_name
   song = re.sub("^\d+\s*-+\s*","", song)	#remove numbers, spaces and dash "10 - Song" ...
   song = re.sub("\s*-\s*","_", song)	#substitute remaining dashes with underscore
   song = re.sub("The %s" %band_name, band_name, song) #remove "The" from the beginning of band names
   song = re.sub("MP3","mp3",song)	#un-capitalize
   song = re.sub("_\s*", "_", song)
   match = re.search("%s_" % band_name,song)  
   
   #sometimes the band name in the database or the newly added already appears with an underscore. we don't want that.
   #so take given band name, add underscores and search for it in the song name. 
   #if exists - remove it (substitute with "") 
   band_name2 = re.sub("\s","_",band_name)   
   match2 = re.search(band_name2, song)
   if (match2): 	#remove name with underscores, if exists.
      song = re.sub(band_name2,"",song)	
   
   #add the band name if it's not already in the filename 
   if (not match ): 
      song = band_name + "_"+song
   else:
      song = re.sub("%s_" %band_name, "", song)
      song = band_name + "_"+song
   
   song = re.sub("_+","_",song)
   #print "working on : ", song 

   try:
      f = MP3(old_name)
      bitrate = f.info.bitrate / 1000
   except:
      print "Critical Error: Problem with file %s" % old_name
      exit() 
   #song is new	because no bit_rate was recorded for it so far in our file-passage (but it still might be there) 
   if (bitrate_dict.get(song)==None):
      bitrate_dict[song] = bitrate
      #name_dict[song] = old_name
      count_dict[song] = 1
      #the name is changed, need to move (rename) file
      #print song,"\n",old_name
      if (song != old_name):
         try:
            if (os.path.exists(song)):
			   #need to decide now between pre-existing version and this one before we rename anything.
			   #compare_songs(br_curr, br_exists, size_curr, size_exists)
               exists_br = MP3(song).info.bitrate/1000
               comp = compare_songs(bitrate, exists_br , os.path.getsize(old_name), os.path.getsize(song))
               #if current song in this loop is better than pre-existing, rename
               if (comp == 1):	
                  print "Deleting pre-existing song in a correct naming system and replacing with a better version!\n"
                  os.rename(old_name, song)
               elif (comp == 2):
                  print "File already existed. Pre-existing version is better: ", exists_br, bitrate
                  bitrate_dict[song] = exists_br			
                  os.remove(old_name) 
            else:
               print "No pre-existing version, renaming..."
               os.rename(old_name, song)  		   
         except:
            print "Warning: failed converting %s to %s" %(old_name, song) 
	     
   #there's already a version of this song in this directory which we encountered and  name-corrected
   else: 
      count_dict[song] += 1
      comp=compare_songs(bitrate, bitrate_dict[song], os.path.getsize(old_name), os.path.getsize(song) )
      if ( comp == 1):
         print old_name," is better: ",bitrate, ">",bitrate_dict[song] 
         bitrate_dict[song] = bitrate
         os.rename(old_name, song)
         print old_name, bitrate, "\t==>\t" , song
		 #name_dict[song] = old_name
      elif ( comp == 2 ):
         #get rid of current fle only if it's smaller
         print song, " previous encountered version is better: ",bitrate_dict[song], ">", bitrate
         os.remove(old_name)
   
for (k,v) in count_dict.items():
   print >>fp, k #str(v)+" : "+k
fp.close()
os.chdir("../")
