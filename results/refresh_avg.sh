#!/bin/bash
echo > time.avg
echo > memory.avg
echo > points.avg
echo > time_point_dist.avg
echo > rate.avg
for f in $(cat conf_avg/f.conf)
do
	for s in $(cat conf_avg/s.conf)
	do
		for t in $(cat conf_avg/t.conf)
		do
			for theta in $(cat conf_avg/theta.conf)
			do
				for l in $(cat conf_avg/l.conf)
				do
					#AVG Time
					sum=0
					nb_tests=0
					for time in $(cat time.all | grep "$f $s $t $theta $l " | cut -d' ' -f6)
					do 
						sum=$(($sum+$time))
						nb_tests=$(($nb_tests+1)) 
					done
					if [ $sum -gt 0 ]
					then
						avg=$(($sum / $nb_tests))
						time_avg=$avg
						echo "$f $s $t $theta $l :$avg: ($nb_tests tests)" >> time.avg
					fi
					
					#AVG Memory
					sum=0
					nb_tests=0
					for mem in $(cat memory.all | grep "$f $s $t $theta $l " | cut -d' ' -f6)
					do 
						sum=$(($sum+$mem))
						nb_tests=$(($nb_tests+1)) 
					done
					if [ $sum -gt 0 ]
					then
						avg=$(($sum / $nb_tests))
						echo "$f $s $t $theta $l :$avg: ($nb_tests tests)" >> memory.avg
					fi
					
					#AVG Points
					sum=0
					nb_tests=0
					for points in $(cat points.all | grep "$f $s $t $theta $l " | cut -d' ' -f6)
					do 
						sum=$(($sum+$points))
						nb_tests=$(($nb_tests+1)) 
					done
					if [ $sum -gt 0 ]
					then
						avg=$(($sum / $nb_tests))
						points_avg=$avg
						echo "$f $s $t $theta $l :$avg: ($nb_tests tests)" >> points.avg
					fi
					
					#AVG Time for one distinguished point
					if [ $sum -gt 0 ]
					then
						time_point=$(($time_avg / $points_avg))
						echo "$f $s $t $theta $l :$time_point: ($nb_tests tests)" >> time_point_dist.avg
					fi
					
					#AVG Rate
					sum=0
					nb_tests=0
					for rate_str in $(cat rate.all | grep "$f $s $t $theta $l " | cut -d' ' -f6)
					do 
						rate=$(echo $rate_str | cut -d'.' -f1)$(echo $rate_str | cut -d'.' -f2)
						sum=$(($sum+$rate))
						nb_tests=$(($nb_tests+1)) 
					done
					if [ $sum -gt 0 ]
					then
						avg=$(($sum / $nb_tests))
						avg_l=$(($avg / 100))
						avg_r=$(($avg % 100))
					fi
                    sum_=0
					nb_tests_=0
                    for rate_str in $(cat rate.all | grep "$f $s $t $theta $l" | cut -d' ' -f7 | cut -d'(' -f2 | cut -d')' -f1)
					do 
                        rate=$(echo $rate_str | cut -d'.' -f1)$(echo $rate_str | cut -d'.' -f2)
                        sum_=$(($sum_+$rate))
                        nb_tests_=$(($nb_tests_+1))
					done
					if [ $sum_ -gt 0 ]
					then
						avg=$(($sum_ / $nb_tests_))
						avg_l_=$(($avg / 100))
						avg_r_=$(($avg % 100))
                        echo "$f $s $t $theta $l :$avg_l.$avg_r ($avg_l_.$avg_r_): ($nb_tests tests)" >> rate.avg
					else
                        if [ $sum -gt 0 ]
				        then
                            echo "$f $s $t $theta $l :$avg_l.$avg_r: ($nb_tests tests)" >> rate.avg
                        fi
                    fi
                    
				done
			done
		done
	done
done
