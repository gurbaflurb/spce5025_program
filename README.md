SPCE 5025 - Homework 9

# Problem 1
Propagate a GEO orbit for 90,000 seconds
- Use a step size of 300, for 300 steps
- Compute ECEF position at each point
- Compute Geodetic latitude, longitude, and altitude at each point
- Plot Latitude Vs. Longitude
![Local Image](./homework9_p1.png)


# Problem 2
Given an Space Vehicle vector, and epoch, and ground site location LION, compute the following:
- Compute Greenwich Hour Angle at vector epoch
- Compute ECEF satellite vector
- Compute ECEF ground site vector
- Compute relative position between SV and site
- Compute ECEF->Topocentric transformation
- Compute relative position in Topocentric frame
- Compute azimuth and elevation of satellite

# Problem 3
For same site and Problem 2 TOD vector
- Compute TOD ground site vector
- Compute range using instantaneous range method
- Using light-time algorithm, compute time to traverse SV->RCVR leg and add /2 to get one-way range
- Let transponder delay time be 1e-6 sec
- Compute difference between instantaneous and light-time range
