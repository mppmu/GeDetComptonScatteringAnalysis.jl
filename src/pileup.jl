baseline_slope(wf) = mean(view(wf, 1400:1500)) - mean(view(wf, 1:100))

function is_not_pileup(wf, cut)
    try
        wf = rm_baseline(wf)
        -cut < mean(view(wf, 1800:2000)) < cut && count_peaks(wf) ≤ 1 
    catch
        false
    end 
end

function is_not_peak_pileup(wf)
    wf = rm_baseline(wf)
    try
        count_peaks(wf) ≤ 1 
    catch
        false
    end 
end

function count_peaks(wf)
    xings_i = []
    xings = 0
    twf = maw(trap_filter(wf,18,250,55.0),100)
    thold = mean(view(twf,1:200)) + 3e6
    for i in 1:length(twf)-1
        if twf[i] ≤ thold && twf[i+1] ≥ thold
            push!(xings_i, i)
        end
    end
    for i in 1:length(xings_i)-1
        if xings_i[i+1] - xings_i[i] > 100
            xings += 1
        end
    end
    xings +=1 #final xing
end

function rm_baseline(wf)
    av = sum(view(wf,1:200))/200
    return wf .- av
end
