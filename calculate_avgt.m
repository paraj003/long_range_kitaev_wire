function ms=calculate_avgt(Ns, Bs, J,narr, x, w)
    % plot <G_per>_t = the time average of G_per(t) vs B.
    %Plot avg for different times given by x(narr)

    ms = zeros(length(Ns), length(Bs),length(narr));
    for j=1:length(Ns)
        m_even = create_m(Ns(j), 0, J, 1);
        m_odd = create_m(Ns(j), 0, J, -1);
        fprintf("N = %d\n", Ns(j));
        for i=1:length(Bs)
            fprintf("Bs = %d\n", Bs(i));
            % Update m with new B value. Don't recreate m. Waste of time.
            m_even = update_m(m_even, Bs(i));
            m_odd = update_m(m_odd, Bs(i));
            ms(j, i,:) = avg_G_per_turbo(m_even, m_odd,narr, x, w);
        end
        mdata=ms(j,:,:);
        N=Ns(j);
        filename=sprintf('data/Avg_mag_N_%G_t_%G:%G:%G_Bs%G:%G:%G.mat',Ns(j),x(1),x(2)-x(1),x(end),Bs(1),(Bs(end)-Bs(1))/(length(Bs)-1),Bs(end))
        save(filename,'x','w','N','Bs','J','narr','mdata')
    end
end

