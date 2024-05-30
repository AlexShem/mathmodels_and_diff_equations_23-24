function Su = S_smoothing(u)
Su = u - (circshift(u, -1) + circshift(u, 1))/2;
end
