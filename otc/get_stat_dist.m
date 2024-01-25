%https://www.mathworks.com/help/econ/dtmc.asymptotics.html#d126e153249
function p = get_stat_dist(P)

mc = dtmc(P);
p = asymptotics(mc);

end
