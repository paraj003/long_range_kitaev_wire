function m = update_m(m, newB)
    % Rather than create a new m everytime, just update what has changed.
    % in this case only the B field.
    for k=1:2:length(m)
        m(k, k) = newB / 2.0;
        m(k+1, k+1) = -newB / 2.0;
    end
end


