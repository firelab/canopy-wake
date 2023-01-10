

The goal here is to try to create a 50 ft x 50 ft x 30 ft (15.2400 m x 15.2400 m x 9.1440 m) cube .stl object, placed with its bottom on the ground in z, with its center straight down the middle of the domain in y, with its front end at a chosen distance D (examples: 4H, 8H, 12H, 16H, 20H) downwind of where a vegetation fence of length H is placed at 10H downwind in x, using snappy hex mesh. The cube is to represent a full physical building, just the geometry is needed for creating the .stl object in freeCAD, boundary conditions for the building set in the /0/ time directory files for snappy hex mesh to work, and snappy hex mesh settings files set in /system/ for snappy hex mesh to work. The size of the building in z comes from assuming ~12 ft per story and assuming the building is ~2 stories tall, which comes out to 24 ft but 30 ft didn't seem too bad and looked easier to work with. The 50 ft x 50 ft size of the building came from estimating that 50 ft x 50 ft to 100 ft x 100 ft would come out to be 2,500 ft^2 to 10,000 ft^2 which seemed to be similar to the dimensions of some small to medium sized houses (googling "average size of a house Los Angeles" gives 2010-2016 New-Home Median size of 1,800 ft^2, Full Housing Stock Median of 1,488 ft^2, values seem to get up to or around 2,500 ft^2).

I took my starting point for this cube and snappy hex mesh stuff, from /1stEmberProjectIdealizedTestCases/1stCases_meshingTechniques/freeCadSnappy_FlowAroundObjects/. The original case was a replication of the example given in this video: https://www.youtube.com/watch?v=KZWHoTpqAOE&list=PLWVpkFXl_Olay-s3kkSzZeRjo_nn7BgDd&index=3&ab_channel=Engineeringtalks, which went over creating .stl objects in freeCAD, followed by making them ascii using meshLab, and finally setting up the surfaceFeatureExtractDict, snappyHexMeshDict, and meshQualityDict files followed by running snappyHexMesh to generate and view the mesh in paraview. It doesn't go over what to do to set the regular files and BCs, but I was able to get the starting BC conditions for the first cases by pausing the video to capture what values were used for BCs in the video.

Anyhow, notice that in the video they have it rotated so when they modify the z position, that is effectively our y position for these cases. It also turns out that you need to be careful to place the objects with the x position at the left side not the center of the desired object, the y position at the center not the bottom side of the desired object. And the positions have to match the desired positions within the domain, so if you change the domain size, the positions of the .stl file generation in freeCAD have to be updated as well.


Phew, it looks like the rest of the snappyHexMesh files didn't really need edited, other than to throw out the extra cylinder stuff. I didn't end up needing to do anything but get the cube object created into the correct position and at the correct size, modify the time directory files to have the BCs for the cube, and then just run it. It did take a while, but it seems to be working fine and it seems to be a decent looking mesh result :). See the full video if need to study more how to improve the snappyHexMesh settings.

Hmm, don't forget to add the patchnames to each variable file in the first time directory. Luckily the ones already figured out for the video case already look good enough, at least for use as a first pass. Looks like you also have to copy and paste the set of files from the 0 directory to the 1 and 2 directory created by snappyHexMesh before running setFields.

BIG NOTE: I had to be careful of units in freeCAD after all, it turns out that you HAVE TO ENTER THE UNITS AS M VALUES, BUT AS MM UNITS. So it takes your m values, and says that they are mm values, but the output comes out to be m values. So if you try to use mm values, the output becomes micrometer values!!!. I was also having trouble when the values are so big that there aren't enough decimal spaces, I was worried that it was cutting off my building locations if the values were over 99.9999 with 4 decimal places as it would convert the units from mm to m to use fewer decimal places in the display, but just letting it do its thing to the number conversion but still using 4 decimal places and the units in m still seemed to give the proper output. There may be a way to mess with this by going to edit->preferences->units, where you can also set decimals, but the little bit I did to try to mess with it just made things worse, (I tried mks at 4 decimals, which confused it worse), so don't trust it until you've thoroughly tested it if you are going to try to mess with that stuff, make sure it behaves consistently when opening and closing the program multiple times with the setting changes and resetting the computer (yes it is that bad, it would act nice till a restart of the program or a reset of the PC, then it would act all screwed up and I had to reinstall to get to what I at least had that was defined behavior).



For this case, I made a new freeCAD part for a cube, sized 50 ft x 50 ft x 30 ft (15.2400 m x 15.2400 m x 9.1440 m). The vegetation fence is 1H x 48H x 1H with H = 6 ft (1.8288 m) representing 48 trees of size L x W x H with L == W == H. So the fence is 6 ft x 288 ft x 6 ft (1.8288 m x 87.7824 m x 1.8288 m).
Now, normally the fence is placed 10H downwind in the domain in x, and I use the fence to figure out the position of the building. But in this case, the domain had to be scaled to the building instead of to the fence as the fence is much smaller than the building, and the building required more distance upwind and downwind than the fence. So in some ways this makes it a heck of a lot easier to place the building.
Anyhow, in x, the building is to be placed at 6Hbdg, so at 6*30 ft = 180 ft (54.8640 m). The y position is a bit trickier, for some odd reason the x needs placed with the left side of the building as the x position, but the y needs placed off to the side of the center of the building. So I need to figure out yMid, then subtract half the building from that to get the y position of the building. So, the y domain size is 6Hbdg + 50 ft (the building) + 6Hbdg = 12Hbdg + 50 ft = 12*30ft + 50 ft = 360 ft + 50 ft = 410 ft (124.968 m). So yMid is 410 ft / 2 = 205 ft (62.484 m). So the y position of the building is 205 ft - 50 ft / 2 = 205 ft - 25 ft = 180 ft (54.8640 m). Oh, this is 6Hbdg, duh. Ah well, the important thing is getting the answer and verifying it. The z position of the building is 0.








To setup the domain with house/building snappy hex mesh for simulation runs, do the following:
  blockMesh 2>&1 | tee blockMesh.log
  checkMesh 2>&1 | tee checkMesh0.log
  surfaceFeatures 2>&1 | tee surfaceFeatures.log
  snappyHexMesh 2>&1 | tee snappyHexMesh.log
  checkMesh 2>&1 | tee checkMesh1.log
  
  Then copy and paste the missing /0/ time directory files to the newly created time directory files (probably /1/ and /2/), and now the mesh should be viewable with paraview/paraFoam before doing a simulation run. Note, sometimes you can copy and paste the whole block of commands and they'll all just run one after the other, showing what they are putting into the log files as they run, putting the runs into the separate log files.















