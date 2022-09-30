

I originally thought that I wouldn't be able to run any of the lagrangian submodels without writing a new lagrangian particle solver that calls them as a new cloud or parcel. Luckily it seems that was not required after all. I then separated out the libraries and worked to get them to compile, then put them all together just to keep a more organized method, then left only the ones that I was actually changing. Trying to get the compile methods right for the dispersion models, cloud function objects, and injection models was extra confusing, but I was able to get it to work, and I was able to test it and see that you can use the old and the new functions/models interchangeably so they are being added to the existing lists of models correctly :). I'm glad it seems to have worked as the guides said to just copy and paste the ENTIRE library and edit the small parts that you want to change, this method is much cleaner.

I originally started out with just a copy of the existing code, with slight name changes, with the goal to just try to compile the dang things. Then I used the original copies as examples for writing my own, making it so this folder only includes functions/objects/models/submodels that are changed from the original copies.


To make this run, just do wmake in the top folder.



