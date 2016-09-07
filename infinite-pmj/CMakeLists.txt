set(PROGS naive_mc order_n)

foreach(name ${PROGS})
  set(target_name ${name})
  add_executable(${target_name} ${name}.cpp)
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${name})
  add_alps_test(${name})
endforeach(name)
