import pandas as pd
import numpy as np
import plotly.express as px
from datetime import datetime
import os

def get_components_of_runtime(table, name="unnamed"):
    sub_tbl = table[["Eviction Time",
                     "Baseline minor PF Time",
                     "Extra Minor PF Time",
                     "Major PF Time",
                     "Baseline User Time",
                     "Extra User Time",
                                     ]] / 1e6
    sub_tbl["Experiment Name"] = table["Experiment Name"]
    fig = px.area(sub_tbl, title='Components of runtime(%s)'%name,
                  color_discrete_sequence=['#636efa', '#ef553b',  '#9e1700','#00cc96', '#ab63fa', '#3c0c73'],
                  animation_frame="Experiment Name")
    fig.update_layout(
        xaxis_title="Ratio",
        yaxis_title="Time(seconds)",
    )
 #   fig.add_trace(px.line(table["Measured(wallclock) runtime"]).data[0])
  #  fig.add_trace(px.line(table["sys+usr"] / 1e6).data[0])

    def anno(text, posx = 1.1, posy=0.32):
        dy = -0.04
        if anno.counter > 0:
            posx += 0.15
        fig.add_annotation(text=text,
              xref="paper", yref="paper",
              x=posx, y=posy + dy * anno.counter, showarrow=False)
        anno.counter+= 1
    anno.counter = 0

    #anno("Workload constants:")
    #anno("Baseline System Time(s): %.2f" % (table["Baseline System Time"].values[0]/1e6))
    #anno("Baseline App Time(s): %.2f" % (table["Baseline App Time(us)"].values[0] / 1e6))
    #anno("Baseline Minor PF Time(us): %.2f" % table["Baseline Single Minor PF Time(us)"].values[0])

    return fig

def get_experiment_data(EXPERIMENT_TYPES, experiment_name, experiment_dir):

    # get experiment data
    get_table = lambda experiment_type, table: pd.read_csv("%s/%s/%s/%s_results.csv" % (experiment_dir, experiment_name, experiment_type, table))
    all_data = pd.DataFrame()
    for exp_type in EXPERIMENT_TYPES:
        cgroup = get_table(exp_type, "cgroup").set_index("RATIO")
        ftrace = get_table(exp_type, "ftrace").set_index("RATIO")

        # in multithreaded apps info is collected per cpu, so let's average it
        # N.B. todo:: does not work well for everything. would be good to add up number of
        # faults,
        ftrace = ftrace.groupby(["RATIO"]).sum()
        time_and_swap = get_table(exp_type, "time_and_swap").set_index("RATIO")
        experiment_final = ftrace.join(cgroup).join(time_and_swap)
        experiment_final["EXP"] = exp_type
        all_data = all_data.append(experiment_final)

    """
    Column Legend
    SWAPIN_* : comes from ftrace, tracks calls to swapin_readahead function and closely measures # of major page faults
    EVICT_*  : comes from ftrace, tracks calls to try_to_free_mem_cgroup_pages. is not used ......
    NUM_FAULTS,NUM_MAJOR_FAULTS: comes from cgroup memory.stat, counts major+minor fault counts
    USER, SYSTEM, WALLCLOCK: from /usr/bin/time
    PAGES_EVICTED,PAGES_SWAPPED_IN : comes from fastswap NIC counters
    """
    return all_data


def augment_tables(tbl, filter_raw=True):

    def per_exp_baseline(exp, measurement_col):
        # for each experiment, val is value of `measurement` under 100% local memory
        val = tbl[(tbl.index == 100) & (tbl["EXP"] == exp)][measurement_col].values[0]
        return val

    tbl["Single Minor PF time(us)"]       = tbl["PAGE_FAULT_TIME"] / tbl["PAGE_FAULT_HIT"]

    for exp in tbl["EXP"].unique():

        # all uppercase names are filtered in the end. Also, "Experiment Name" is more
        # readable for plot sliders
        tbl.loc[tbl["EXP"] == exp, "Experiment Name"] = exp

        # set appropriate baseline column, filter rows by exp type(tbl["EXP"] == exp)
        tbl.loc[tbl["EXP"] == exp, "Baseline System Time"] = per_exp_baseline(exp, "SYSTEM") * 1e6
        tbl.loc[tbl["EXP"] == exp, "Baseline User Time"] =   per_exp_baseline(exp, "USER") * 1e6
        tbl.loc[tbl["EXP"] == exp, "Baseline Single Minor PF Time(us)"] =   per_exp_baseline(exp, "Single Minor PF time(us)")


    tbl["Page Faults"]             = tbl["PAGE_FAULT_HIT"].fillna(0)
    tbl["PF Time(us)"]             = tbl["PAGE_FAULT_TIME"].fillna(0)
    # can also be "NUM_MAJOR_FAULTS"
    tbl["Major Page Faults"]       = tbl["SWAPIN_HIT"].fillna(0)
    tbl["Minor Page Faults"]       = tbl["NUM_FAULTS"] - tbl["NUM_MAJOR_FAULTS"]
    tbl["Evictions"]               = tbl["PAGES_EVICTED"].fillna(0)

    # deprecate this..does not consider fastswap offload..
    #tbl["Eviction Time(us)"]      = tbl["Evictions"] * SYSTEM_CONSTANTS["SYSTEM_TIME_PER_EVICTION"]
    tbl["Eviction Time"]           = tbl['EVICT_TIME'].fillna(0)

    tbl["Major PF Time"]       = tbl["SWAPIN_TIME"].fillna(0)

    # todo:: sth wrong with Sync time since:
   #    it is called by do_page_fault and its time should be included in minor fault time
   #    it is not present in linux_prefetching baseline
   #    ----- so, extra minor PF time should be AT LEAST sync time
   # conclusion does not hold since at times SYNC time is 0.6 sec and extra minor fault time is 0.5
   # tbl["Tape Sync Time(us)"]      = tbl["SYNC_TIME"].fillna(0)
    tbl["Baseline minor PF Time"]  = tbl["Minor Page Faults"] * tbl["Baseline Single Minor PF Time(us)"]

    tbl["Extra Minor PF Time"] = (tbl["PF Time(us)"] - tbl["Major PF Time"] -  tbl["Baseline minor PF Time"]).clip(0)

    tbl["Extra User Time"]         = tbl["USER"] * 1e6 - tbl["Baseline User Time"]

    tbl["System Overhead"]         = tbl["Baseline minor PF Time"] + \
                                       tbl["Extra Minor PF Time"] + \
                                       tbl["Eviction Time"] + \
                                       tbl["Major PF Time"]

    tbl["Total System Time"]       = tbl["System Overhead"] + \
                                       tbl["Baseline System Time"]

    def to_seconds(a):
        pt = datetime.strptime(a,'%M:%S.%f')
        total_seconds = pt.microsecond * 1e-6 + pt.second + pt.minute*60 + pt.hour*3600
        return total_seconds
    tbl["Measured(wallclock) runtime"] = tbl["WALLCLOCK"].map(to_seconds)

    tbl["Runtime"]                 = tbl["Measured(wallclock) runtime"]
    tbl["sys+usr"]                 = tbl["USER"] * 1e6 + tbl["SYSTEM"] * 1e6
    tbl["Runtime w/o Evictions"]   = tbl["Runtime"] - tbl["Eviction Time"]/1e6

    degr = lambda exp, baseline, c: tbl.loc[tbl["EXP"] == exp, c] / baseline

    for exp in tbl["EXP"].unique():
        baseline = per_exp_baseline(exp, "Runtime")
        tbl.loc[tbl["EXP"] == exp, "Degradation(%)"] = degr(exp, baseline, "Runtime") * 100

        baseline = per_exp_baseline(exp, "Runtime w/o Evictions")
        tbl.loc[tbl["EXP"] == exp, "Degradation w/o Evictions(%)"] = degr(exp, baseline, "Runtime w/o Evictions") * 100

        tbl.loc[tbl["EXP"] == exp, "Baseline App Time(us)"] = per_exp_baseline(exp, "Measured(wallclock) runtime")

    if filter_raw:
        raw_cols = [c for c in tbl.columns.values if c.upper() == c]
        tbl.drop(columns=raw_cols, inplace=True)

    return tbl

def take_column_named(column_name, data):
    res = pd.DataFrame()

    for name,df in data.items():

        ## TODO::::: VERY VERY HACKY.. assumes all tables have data appropriately sorted
        ## sanity check later with more data, if these plots become crucial
        if "Experiment Name" not in res:
            res["Experiment Name"] = df["Experiment Name"]

        res["(%s)%s" % (name,column_name)] = df[column_name]
    return res

def get_nic_monitor_data(EXPERIMENT_TYPES, workload, experiment_dir):
    all_data = pd.DataFrame()
    for exp in EXPERIMENT_TYPES:
        tbl = pd.DataFrame()
        for f in os.listdir("%s/%s/%s/" % (experiment_dir, workload, exp)):
            if f.startswith("nic_monitor"):
                ratio = int(f.split(".")[1])
                tmp = pd.read_csv("%s/%s/%s/%s" % (experiment_dir, workload, exp, f))
                tmp["RATIO"] = ratio
                tmp["Experiment Name"] = exp
                tbl = tbl.append(tmp)
        all_data = all_data.append(tbl)

        all_data["Time(s)"] = all_data["TIME"] / 1000
        all_data["Xmit(MB)"] = all_data["XMIT"] / (1024 * 1024)
        all_data["Recv(MB)"] = all_data["RECV"] / (1024 * 1024)

    return all_data.sort_values(["RATIO", "TIME"])
