-- .nvim.lua — Zig project local config
-- This file is trusted per-project; load with :set exrc

-- .nvim.lua (project root)
vim.api.nvim_create_autocmd("FileType", {
  pattern = "zig",
  callback = function()
    vim.bo.makeprg = "time zig test src/main_test.zig --summary all"
  end,
})
