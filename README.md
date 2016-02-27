#wx_scraper
##python code to fetch and plot station data from Weather Underground

###Info
http://areed145.github.io/wx_scraper

http://www.adammreeder.com

http://www.flickr.com/people/adamreeder

### The Reason
If you, like me, contribute station data to Weather Underground, you might be frustrated with how difficult it can be to access higher resolution station data over long periods of time. As an owner of a station that uploads exclusively to WU, I have no other way to access the data.

My interest in actually getting the raw(ish) data from my station grew from a desire to create wind rose plots. As a pilot of aircraft and paragliders, winds are just particularly interesting to me. For some strange reason this curiosity seemed worth the effort.

### The Solution
WU allows downloads of raw(ish) data for single days so the strategy here is to hit their API iteratively over the date range, and store that data locally. After the initial fetch, the program updates only what days aren't stored locally plus the last 7 days of data to make sure it's fresh.

The plots are somewhat old-school in feel, that's by design. The critical ones to me were obviously the wind rose plots, and I'm a fan of the combo plots as well.

There are a few calculations incorporated:
 - Minimum Cloud Base (ft) = ((TempF - DewpointF) / 4.4) * 1000 + Elevation
    - in theory this is the lowest clouds can exist assuming normal adiabatic lapse rate over the station
 - dTdt (degF/hr) = (T1 - T2)/(t1 - t2)
    - rate of temp change; plotted vs solar radiation (W/m^2) and used in combo plots
 - dPdt (inHg/hr) = (P1 - P2)/(t1 - t2)
    - rate of pressure change; not used at the moment

Let me know of any improvement ideas! I enjoy working on this and plan to keep it updated.

Thanks!
