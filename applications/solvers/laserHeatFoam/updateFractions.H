// Consolidated fraction records the maximum liquid fraction ever reached.
rc = max(rc, fL);

// Powder, melt, and re-solidified fractions.
rp = scalar(1) - rc;
rm = fL;
rs = rc - fL;

// Binary history field for cells that have reached the liquidus.
meltHistory = max(meltHistory, pos(T - Tl));
