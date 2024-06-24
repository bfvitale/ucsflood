#!/opt/python-3.11.5/bin/python3

import regionlib

all_scenarios = []

for year in range(2020, 2100 + 1, 10):
    for severity in [ regionlib.Severity.INT_LOW,
                      regionlib.Severity.INT,
                      regionlib.Severity.HIGH ]:
        for freq in [ 2, 12, 26 ]:
            exparms = regionlib.ExParms(year, severity, freq)
            all_scenarios.append(regionlib.exparms_suffix(exparms))


priority_by_scenario = { scenario: int(priority)
                         for priority, scenario in [
                                 line.strip().split()
                                 for line in open('priority-ucs')] }
for scenario in all_scenarios:
    if scenario not in priority_by_scenario:
        priority_by_scenario[scenario] = len(priority_by_scenario)

items_sorted = sorted(priority_by_scenario.items(), key=lambda kv: kv[1])
for scenario, priority in items_sorted:
    for region in [
            regionlib.Region.GU, regionlib.Region.PRVI,
            regionlib.Region.WEST, regionlib.Region.EAST]:
        print(scenario, region)
