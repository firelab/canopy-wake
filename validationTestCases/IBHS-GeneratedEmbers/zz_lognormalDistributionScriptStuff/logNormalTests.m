

mean = -3.2811;
variance = 0.94289;

nPars = 50016;

%%% lognormal distributions

pd = makedist('Lognormal','mu',mean,'sigma',variance);

% standard full lognormal distribution plot for this mean and variance
r = random(pd,nPars,1);
h = histogram(r);
xlim([0,0.5]);


% what I thought was the inputs, my truncation limits
%%% this doesn't look like the excel version at all
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);
y = pdf(pd,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal 0.001 to 0.03")
r = random(pd,nPars,1);
h = histogram(r);
xlim([minX,maxX]);

% what resulted as the limits for the truncated lognormal distribution
%%% this doesn't look like the excel version at all
minX = 0.002675;
maxX = 0.029999;
nx = 1000;
x = linspace(minX,maxX,nx);
y = pdf(pd,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal 0.002675 to 0.029999")
r = random(pd,nPars,1);
h = histogram(r);
xlim([minX,maxX]);

% what resulted as the limits for the nontruncated lognormal distribution
%%% so far this is the only one close to matching the excel plots, the truncation
%%% must do weird stuff to it. There must be an error that is letting the
%%% normal distribution act like it is greater than 0 cause here I'm not
%%% getting greater than 0. Using a positive mean causes the lognormal
%%% distribution to go crazy big, so the negative number appears to be
%%% correct, for a normal distribution it shifts way too much to the right
%%% for the positive mean to be correct, must be some other quirk that lets
%%% a negative mean of this value work for the normal distribution
minX = 0.002005;
maxX = 0.513436;
nx = 1000;
x = linspace(minX,maxX,nx);
y = pdf(pd,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal 0.002005 to 0.513436")
r = random(pd,nPars,1);
h = histogram(r);
xlim([minX,maxX]);



%%% normal distributions

% standard full normal distribution plot for this mean and variance
%%% interesting, when I made this a different way, I got values, even if
%%% just a few, past the 0 value, this is not the case here
r = random('Normal',mean,variance,nPars,1);
h = histogram(r);
%xlim([0,0.5]);


% what I thought was the inputs, my truncation limits
%%% this doesn't look like the excel plot at all
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);
y = normpdf(x,mean,variance);
plot(x,y);
xlabel('x');
ylabel('p');
title("normal 0.001 to 0.03")
r = random('Normal',mean,variance,nPars,1);
h = histogram(r);
xlim([minX,maxX]);

% what resulted as the limits for the truncated normal distribution
%%% this doesn't look like the excel plots at all
minX = 0.001;
maxX = 0.026683;
nx = 1000;
x = linspace(minX,maxX,nx);
y = normpdf(x,mean,variance);
plot(x,y);
xlabel('x');
ylabel('p');
title("normal 0.001 to 0.026683")
r = random('Normal',mean,variance,nPars,1);
h = histogram(r);
xlim([minX,maxX]);

% what resulted as the limits for the nontruncated normal distribution
%%% this doesn't look like the excel plot at all
minX = 0.00000001;
maxX = 1.39801695;
nx = 1000;
x = linspace(minX,maxX,nx);
y = normpdf(x,mean,variance);
plot(x,y);
xlabel('x');
ylabel('p');
title("normal 0.00000001 to 1.39801695")
r = random('Normal',mean,variance,nPars,1);
h = histogram(r);
xlim([minX,maxX]);
xlim([0,0.05]);


%%%%%%%%%% redoing the above best versions, but now trying to truncate the
%%%%%%%%%% distributions to see if that improves the sampling. These look 
%%%%%%%%%% a LOT closer to what I saw in the excel plots, and it is a
%%%%%%%%%% better method for generating normal distributions too

%%% truncated lognormal

pd = makedist('Lognormal','mu',mean,'sigma',variance);

% standard full lognormal distribution plot for this mean and variance
r = random(pd,nPars,1);
h = histogram(r);


% what I thought was the inputs, my truncation limits
%%% hrm, this seems to have elements that are SIMILAR to the excel plot,
%%% especially if you flatten this plot by shrinking the plot window in the
%%% y direction. But it doesn't seem to have the same slope of increasing
%%% value, it seems to happen more smoothly rather than suddenly compared
%%% to the excel plot. It also decays more after the middle, unlike the
%%% excel plot.
%%% oh interesting, if you use nBins 48 for the histogram, the same number
%%% of bins for the excel histogram, then you get the weird clip off at the
%%% end. Make the window crazy small in y for this, and the steepness of
%%% the plot looks a lot more similar to the excel version. The steepness
%%% looks better for this version than when I truncate it by the resulting
%%% limits rather than these limits.
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal 0.001 to 0.03")

r = random(t,nPars,1);
%h = histogram(r);
h = histogram(r,48);


% what resulted as the limits for the truncated lognormal distribution
%%% hrm, this seems to have elements that are SIMILAR to the excel plot,
%%% especially if you flatten this plot by shrinking the plot window in the
%%% y direction. But it doesn't seem to have the same slope of increasing
%%% value, it seems to happen more smoothly rather than suddenly compared
%%% to the excel plot. It also decays more after the middle, unlike the
%%% excel plot.
%%% oh interesting, if you use nBins 48 for the histogram, the same number
%%% of bins for the excel histogram, then you get the weird clip off at the
%%% end. Make the window crazy small in y for this, and the steepness of
%%% the plot looks a lot more similar to the excel version. The steepness
%%% doesn't look as good for this as it does for the one above where the
%%% limits are the original limits.
minX = 0.002675;
maxX = 0.029999;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal 0.002675 to 0.029999")

r = random(t,nPars,1);
h = histogram(r);
%h = histogram(r,48);

% what resulted as the limits for the nontruncated lognormal distribution
%%% MAN this looks so close to the same as the excel plot. Especially if
%%% you shrink the window in the y direction. Especially if you use the
%%% same number of histogram bins (50), it even gets the upper count close
%%% to the same, though the first initial one is a bit lower in this than
%%% the excel histogram plot when doing it that way.
%%% too bad this isn't grabbing particles at a smaller value
minX = 0.002005;
maxX = 0.513436;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal 0.002005 to 0.513436")

r = random(t,nPars,1);
h = histogram(r);
%h = histogram(r,50);

% trying out with infinite limits rather than the resulted infinite limits just above
%%% Seems to still cut out particles from the range of 0.001 to 0.002, sad
%%% because it seems close enough to what is wanted other than that.
minX = 0.00000001;
maxX = 9999999999;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal inf limits")

r = random(t,nPars,1);
h = histogram(r);
xlim([0,0.5]);


%%% truncated normal

pd = makedist('Normal','mu',mean,'sigma',variance);

% standard full normal distribution plot for this mean and variance
%%% interesting, this method appears to be better than the other one for
%%% generating the normal distribution. There aren't many, but in this one
%%% there are just a few values past 0, which was not true with the past
%%% method
r = random(pd,nPars,1);
h = histogram(r);


% what resulted as the limits for the truncated normal distribution
%%% this one doesn't look like the excel plot of it. Well, if you clip off
%%% the edges from both plots, it seems like it would look similar. So it
%%% could just be some annoying difference in how the bins are calculated
%%% for the separate software, looking at the bin limits of excel, they are
%%% off by half a cell from the bin limits of this plot.
minX = 0.001;
maxX = 0.026683;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("normal 0.001 to 0.026683")

r = random(t,nPars,1);
h = histogram(r);

% what resulted as the limits for the nontruncated normal distribution
%%% WOW this one looks quite similar to the excel one even though the just
%%% above doesn't look the same. Using 50 bins for the histogram does
%%% increase the nparticles per bin, but not enough to be the same numbers
%%% as the excel plots though.
%%% having a tough time telling whether my lognormal distribution works
%%% correctly or not, it seems to be working correctly from what I can tell
%%% now, or at least close enough, but it doesn't seem to be able to catch
%%% the extra small particles that are caught by this distribution.
%%% turns out that the problem is that the distribution is still flat at
%%% the range of 0.001 m to 0.03 m, just APPEARS to be lognormal over this
%%% bigger range. Even if it has a similar but not the same slope as
%%% lognormal, there is another drop to the left of the peak in lognormal
%%% that isn't there with these shifted normal distributions.
minX = 0.00000001;
maxX = 1.39801695;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("normal 0.00000001 to 1.39801695")

r = random(t,nPars,1);
h = histogram(r);
%h = histogram(r,50);



%%
%%%%%% try out a truncated normal distribution with a different mean and
%%%%%% variance, see if that looks better

nPars = 50016;

%%% try mean of 0, same variance
mean = 0;
variance = 0.94289;

pd = makedist('Normal','mu',mean,'sigma',variance);

% standard full normal distribution plot for this mean and variance
%%% seems like it may work
r = random(pd,nPars,1);
h = histogram(r);

% looks like this one is like a slightly shifted uniform distribution, so
% not good enough after all.
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("mean 0 variance 0.94289")

r = random(t,nPars,1);
h = histogram(r);


%%% try mean of 0, variance of 0.01
mean = 0;
variance = 0.01;

pd = makedist('Normal','mu',mean,'sigma',variance);

% standard full normal distribution plot for this mean and variance
%%% seems like it may work
r = random(pd,nPars,1);
h = histogram(r);

% wow, this one actually looks relatively close. Turns out I thought it
% didn't before, because I forgot to update the pdf and used the same pdf
% for all the cases. Anyhow, this one has a different slope than the
% lognormal case, and it probably could use more particles on the higher
% end than is given here, but it is surprisingly close
% turns out it is also missing the lower number of particles to the left of
% the curve.
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("mean 0 variance 0.01")

r = random(t,nPars,1);
h = histogram(r);


%%% try mean of -2, same variance as original
mean = -2;
variance = 0.94289;

pd = makedist('Normal','mu',mean,'sigma',variance);

% standard full normal distribution plot for this mean and variance
%%% looks a lot more like the tail end, but the tail is now shifted further
%%% over, so may work
r = random(pd,nPars,1);
h = histogram(r);

% still looks flat, dang
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("mean -2 variance 0.94289")

r = random(t,nPars,1);
h = histogram(r);


%%%% yeah, it appears that it has to be the original, but not truncated 
%%%% limits for it to work, not the original truncated limits. This is 
%%%% because the original truncated limits is too uniform and does not 
%%%% appear to be lognormal at all, but the original not truncated limits 
%%%% appears to be lognormal.
%%%% But this is misleading because the original not truncated limits turns
%%%% out to not look lognormal over the range of 0.001 m to 0.03 m sized
%%%% particles that we are interested in, even though it looks lognormal
%%%% over the full not truncated range. Even though it is lognormal-like
%%%% over the full not truncated range, it is still uniform-like over the
%%%% desired range. Could try to get around this by using the not truncated
%%%% range, but then the desired range within the range is uniform, not
%%%% lognormal, and the not truncated range introduces a ton more big 
%%%% particles that are bigger than 0.03 m. I don't think that we want 
%%%% those particles.


%%% hrm, the part of the normal distribution not truncated that we seem to
%%% like, goes from 0.0001 m to 0.61 m, so if I try a mean of original
%%% value - 0.61 m, what does that look like?

mean = -3.2811-0.61;
variance = 0.94289;

pd = makedist('Normal','mu',mean,'sigma',variance);

% standard full normal distribution plot for this mean and variance
%%% looks a lot more like the tail end, but the tail is now shifted further
%%% over, so may work
r = random(pd,nPars,1);
h = histogram(r);

% still looks like a slowly decreasing uniform distribution, not so flat as
% before, but it is still too flat, dang
minX = 0.00001;
maxX = 0.1;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("mean -3.8911 variance 0.94289")

r = random(t,nPars,1);
h = histogram(r);


%%% doesn't seem to be right either, dang. I'm starting to think that
%%% lognormal behavior is actually some kind of plotting illusion, from
%%% using really long x values but not so big y values or something.


%%% so I looked at the wiki on normal distribution, how different mean and
%%% variance changes the distribution, and it looks like a different mean
%%% just shifts it left or right, a variance of 1 is a standard amount of
%%% steepness, a variance larger than 1 flattens the curve, a variance
%%% smaller than 1 steepens the curve
%%% part of the problem with the above, is that while doing a nontruncated
%%% normal distribution looks like it is lognormal, over the range of 0.001
%%% m to 0.03 m that we care about, it is still mostly a uniform
%%% distribution, I'm trying to change that.
%%% so I'm going to now try a mean of 0, a but varying the variance

% mean of 0, variance of 0.1

mean = 0;
variance = 0.1;

pd = makedist('Normal','mu',mean,'sigma',variance);

% standard full normal distribution plot for this mean and variance
%%% looks a lot steeper, so hopefully this will work
r = random(pd,nPars,1);
h = histogram(r);

% still looks flat, dang
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("mean 0 variance 0.1")

r = random(t,nPars,1);
h = histogram(r);


% mean of 0, variance of 5

mean = 0;
variance = 5;

pd = makedist('Normal','mu',mean,'sigma',variance);

% standard full normal distribution plot for this mean and variance
%%% looks a lot wider, not sure that this will work
r = random(pd,nPars,1);
h = histogram(r);

% still looks flat, dang, actually even more flat
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("mean 0 variance 5")

r = random(t,nPars,1);
h = histogram(r);


%%% looks like the steeper the better, but it needs to be crazily steeper
%%% than what I had before. I had one that seemed close with a variance of
%%% 0.01 before, but what does it look like when using an even smaller
%%% variance?

% mean of 0.001 (had to affect cause the variance made the range too small),
% variance of 0.0001

mean = 0.001;
variance = 0.0001;

pd = makedist('Normal','mu',mean,'sigma',variance);

% standard full normal distribution plot for this mean and variance
%%% looks a lot steeper, so hopefully this will work
r = random(pd,nPars,1);
h = histogram(r);

% still looks flat, dang, actually even more flat
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("mean 0.001 variance 0.0001")

r = random(t,nPars,1);
h = histogram(r);


%%%%%
%%% I had tried an even smaller variance, and it broke because the data was
%%% not on the line. Looks like I need to use the mean 0 and variance 0.01
%%% as the starting point and go from there, that was the best case and was
%%% surprisingly close

% mean of 0 variance of 0.01. This is a repeat of the best case above

mean = 0;
variance = 0.01;

pd = makedist('Normal','mu',mean,'sigma',variance);

% standard full normal distribution plot for this mean and variance
%%% looks a lot steeper, so hopefully this will work
r = random(pd,nPars,1);
h = histogram(r);

% the steepness seems relatively okay, but the range of values needs more
% data on the right side of the distribution. Not sure if this means it
% needs something a bit more steep, or just to shift the whole thing to the
% right a bit. Technically the log normal has the opposite direction slope
% on the left side as well, which wouldn't be caught here without shifting
% the thing a bunch to the right.
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("mean 0 variance 0.01")

r = random(t,nPars,1);
h = histogram(r);


%%% try shifting it to the right a little bit, not trying to capture the
%%% left side of the lognormal yet though, even though this SHOULD capture
%%% that left side a little bit

% mean of 0.003 variance of 0.01.

mean = 0.003;
variance = 0.01;

pd = makedist('Normal','mu',mean,'sigma',variance);

% standard full normal distribution plot for this mean and variance
%%% looks a lot steeper, so hopefully this will work
r = random(pd,nPars,1);
h = histogram(r);

% Looks like it needs shifted even further to the right, it also seems like
% it needs to be steeper on the left and less steep on the right to get the
% lognormal look
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("mean 0.003 variance 0.01")

r = random(t,nPars,1);
h = histogram(r);


%%% try shifting to the right a little bunch more, no change to steepness
%%% yet

% mean of 0.01 variance of 0.01.

mean = 0.01;
variance = 0.01;

pd = makedist('Normal','mu',mean,'sigma',variance);

% standard full normal distribution plot for this mean and variance
%%% looks a lot steeper, so hopefully this will work
r = random(pd,nPars,1);
h = histogram(r);

% Wow! This gets a lot more values for the right side bigger values, but
% also seems to have the dipping off to the left side. It does still seem
% like it needs to be less steep on the right side and more steep on the
% left side to truly match a lognormal profile though.
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("mean 0.01 variance 0.01")

r = random(t,nPars,1);
h = histogram(r);


%%% try the same as before, but now less steep

% mean of 0.01 variance of 0.012 or 0.015.

mean = 0.01;
variance = 0.012;

pd = makedist('Normal','mu',mean,'sigma',variance);

% standard full normal distribution plot for this mean and variance
%%% looks a lot steeper, so hopefully this will work
r = random(pd,nPars,1);
h = histogram(r);

% LOL even a variance of 0.02 was too much, and it became an almost flat
% uniform profile. I can't decide which is better, a variance of 0.012 or
% 0.015. 0.015 seems to get the right steepness for part of the right side,
% loses too much steepness on the left side and the left side of the right
% side, but that may not be a bad thing. They both seem pretty good to me,
% though I'm pretty sure lognormal will always be a better shape.
%%% this is the best set of cases that I've found so far for a normal
%%% distribution. If the lognormal ends up not working as expected, this is
%%% probably what I should use. But if I can get the lognormal to work
%%% correctly, that is what I should use.
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("mean 0.01 variance 0.012")

r = random(t,nPars,1);
h = histogram(r);



%% the true lognormal stuff is being done here

nPars = 50016;

%%% for the lognormal plots, it looks like the mean just shifts it left or
%%% right, where this time the mean is actually the shift from 0 instead of
%%% choosing the middle value, so it is like choosing where to place the
%%% left tail. It looks like the smaller the variance the steeper the
%%% curve, to the point where a variance of 0.25 makes it look like a very
%%% steep normal distribution. Values less than 1 seems to make it go up,
%%% values less than 1 seem to flatten it out more, also shifting the curve
%%% to the left a bunch more as a larger value lets the curve cover a
%%% larger area.
%%% original mean and variance were -3.2811 and 0.94289. Looks like this
%%% means it is shifted almost completely to the left near zero, but at the
%%% same time, it is weird because it is crazy steep even though the
%%% variance is 0.94289. I guess this is because it doesn't seem to let it
%%% shift to the left even though the variance is negative, it just shifts
%%% the middle of the bell closer to zero by a huge amount.
%%% I'm having trouble deciding what needs to be done, the original just
%%% looking at the interval of interest, it seems to cut off a TON of the
%%% tail to the right, and not having enough stuff to the left between 0
%%% and 0.001 m.
%%% I guess this means that I need to either increase the variance or
%%% decrease the mean even more. To retain the curve to the right, does
%%% this mean that I need to do both a little bit at once? Hard to tell.

% original

mean = -3.2811;
variance = 0.94289;

pd = makedist('Lognormal','mu',mean,'sigma',variance);

% standard full lognormal distribution plot for this mean and variance
r = random(pd,nPars,1);
h = histogram(r);
%xlim([0,0.05]);

% over the desired interval, seems to be quite flat for the right hand side
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal mean -3.2811 variance 0.94289")

r = random(t,nPars,1);
h = histogram(r);



% mean of -4, variance of 1

mean = -4;
variance = 1;

pd = makedist('Lognormal','mu',mean,'sigma',variance);

% standard full lognormal distribution plot for this mean and variance
r = random(pd,nPars,1);
h = histogram(r);
%xlim([0,0.05]);

% looks to get it more closely shifted to the left, still could probably go
% to the left a little more, but it also seems to have made the right side
% of it become much more steep. I guess this means needing a bigger
% variance? Not sure, try it with a smaller variance.
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal mean -4 variance 1")

r = random(t,nPars,1);
h = histogram(r);


% mean of -4, variance of 0.8, Try a smaller variance

mean = -4;
variance = 0.8;

pd = makedist('Lognormal','mu',mean,'sigma',variance);

% standard full lognormal distribution plot for this mean and variance
r = random(pd,nPars,1);
h = histogram(r);
%xlim([0,0.05]);

% this looks closer to the original, lost the benefit of stuff going to the
% left as much, but also compared to the original, it is still a heck of a
% lot steeper on the right. I guess I need to try the same but an even
% greater variance.
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal mean -4 variance 0.8")

r = random(t,nPars,1);
h = histogram(r);


% mean of -4, variance of 1.2, Try a greater variance

mean = -4;
variance = 1.2;

pd = makedist('Lognormal','mu',mean,'sigma',variance);

% standard full lognormal distribution plot for this mean and variance
r = random(pd,nPars,1);
h = histogram(r);
%xlim([0,0.05]);

% wow this one actually seems pretty good. Got a lot more on the left side.
% Well, the right side is still too steep, dang. Try a slightly smaller
% mean but greater variance
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal mean -4 variance 1.2")

r = random(t,nPars,1);
h = histogram(r);


%%% try a greater mean and still large variance
% mean of -3.5, variance of 1.2

mean = -3.5;
variance = 1.2;

pd = makedist('Lognormal','mu',mean,'sigma',variance);

% standard full lognormal distribution plot for this mean and variance
r = random(pd,nPars,1);
h = histogram(r);
%xlim([0,0.05]);

% seems even closer. Still enough pars to the left, got steeper on the
% right, though still not enough. Try the same but increase the variance
% even more.
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal mean -3.5 variance 1.2")

r = random(t,nPars,1);
h = histogram(r);


%%% try same mean and still larger variance
% mean of -3.5, variance of 1.4

mean = -3.5;
variance = 1.4;

pd = makedist('Lognormal','mu',mean,'sigma',variance);

% standard full lognormal distribution plot for this mean and variance
r = random(pd,nPars,1);
h = histogram(r);
%xlim([0,0.05]);

% seems to not be getting much unsteeper. was flatter with the last set of
% values. Man this is confusing!!!
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal mean -3.5 variance 1.4")

r = random(t,nPars,1);
h = histogram(r);


%%% try bigger mean, original variance

mean = -3.5;
variance = 0.94289;

pd = makedist('Lognormal','mu',mean,'sigma',variance);

% standard full lognormal distribution plot for this mean and variance
r = random(pd,nPars,1);
h = histogram(r);
%xlim([0,0.05]);

% seems close, need more particles to the left though. Seems to be starting
% to lose the flatness on the right curve though.
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal mean -3.5 variance 0.94289")

r = random(t,nPars,1);
h = histogram(r);



%%% BLARG I'M JUST GOING TO TRY A BUNCH OF STUFF TILL I FIND SOMETHING
%%% THAT LOOKS RIGHT

mean = -3.2811;
variance = 1.2;

pd = makedist('Lognormal','mu',mean,'sigma',variance);

% standard full lognormal distribution plot for this mean and variance
r = random(pd,nPars,1);
h = histogram(r);
%xlim([0,0.05]);

% this seems to have come out to be the closest, ironically keeping the
% original mean, and increasing the variance, but not too far till it went
% crazy, seemed to be what worked.
% seems to be more values to the left, which is good, more of the smaller
% particles. A lot steeper, much more curved down, losing some of the bigger
% particles on the right, but hopefully it isn't too much.
% seems like this is probably the closest thing to what we should use, IF
% the original idea of log profile was even good to begin with.
% this is what we need to use.
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal mean -3.2811 variance 1.2")

r = random(t,nPars,1);
h = histogram(r);


%%% that one was surprisingly close. But when I put it into the code, still
%%% not quite enough particles to the left. Going to keep the variance the
%%% same but now decrease the mean even a little more.

mean = -4;
variance = 1.2;

pd = makedist('Lognormal','mu',mean,'sigma',variance);

% standard full lognormal distribution plot for this mean and variance
r = random(pd,nPars,1);
h = histogram(r);
%xlim([0,0.05]);

% this seems to have come out to be the closest, ironically keeping the
% original mean, and increasing the variance, but not too far till it went
% crazy, seemed to be what worked.
% seems to be more values to the left, which is good, more of the smaller
% particles. A lot steeper, much more curved down, losing some of the bigger
% particles on the right, but hopefully it isn't too much.
% seems like this is probably the closest thing to what we should use, IF
% the original idea of log profile was even good to begin with.
% this is what we need to use.
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal mean -4 variance 1.2")

r = random(t,nPars,1);
h = histogram(r);



%% try fitdist for some of the datasets, see how they compare

mean = -3.2811;
variance = 0.94289;

pd = makedist('Lognormal','mu',mean,'sigma',variance);

% standard full lognormal distribution plot for this mean and variance
r = random(pd,nPars,1);
h = histogram(r);
%xlim([0,0.05]);

% this is another copy of the original
minX = 0.001;
maxX = 0.03;
nx = 1000;
x = linspace(minX,maxX,nx);

t = truncate(pd,minX,maxX);

y = pdf(t,x);
plot(x,y);
xlabel('x');
ylabel('p');
title("lognormal mean -3.2811 variance 0.94289")

r = random(t,nPars,1);
h = histogram(r);


% these commented out lines require someone to put
% diameters = 0;
% into the terminal, then in the workspace, open diameters
% and copy and paste all the diameters into the container
% then that variable can be used here
diameters = 0;


pd_lognormal = fitdist(diameters,'Lognormal')
pd_normal =  fitdist(diameters,'Normal')

% hm, the above doesn't seem to be working. Using it for the different data
% gets me the following values:
% truncated lognormal: mean = -3.96515   [-3.96833, -3.96197], variance = 0.362622   [0.360389, 0.364884]
% truncated normal: mean = 0.0118326   [0.0117621, 0.0119031], variance = 0.00804403   [0.00799449, 0.00809419]
% not truncated lognormal: mean = -3.2834   [-3.28927, -3.27753], variance = 0.669381   [0.665258, 0.673555]
% not truncated normal: mean = 0.123336   [0.12227, 0.124403], variance = 0.121693   [0.120944, 0.122452]
% the only one that gets somewhat close, is the not truncated lognormal.
% But they all look pretty bad. I guess technically it is a set of values
% modified by the min and max desired values, not just a set of values
% going from 0 to 1 modified by the mean and variance.
% guess I'll just have to compare the pdf output to the excel version of
% it. Eww. Or, do a histogram of it and compare that.

h = histogram(diameters);

%%% hrm, I saved these figures, and they looked similar to the excel ones,
%%% but much better. I then went back to where I did the histogram stuff,
%%% and I got a very different set of histograms. Not sure how to tell
%%% whether the lognormal stuff is working correctly or not, it does seem
%%% like there can be artifacts or errors when doing a given set of
%%% distributions, that somehow doesn't happen when replicating the
%%% distribution here in matlab. Grr, so what is the correct thing to do?
%%% I guess I need to make an attempt at a choice of a distribution, and
%%% see if it makes sense when I plot it. At least in here it gave me an
%%% idea of what values to use. Something is still off, something in how it
%%% takes the distribution and turns it into diameters.
%%% does seem like the normal distributions get turned into something weird
%%% looking more than the lognormal distributions. Definitely seems to
%%% match the excel histograms quite well, which is a good sign. Can
%%% probably just keep trying to use the excel stuff.



